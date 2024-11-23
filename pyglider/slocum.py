"""
Routines to convert raw slocum dinkum files to netcdf timeseries.

"""
import bitstring
from datetime import datetime
try:
    import dbdreader
    have_dbdreader = True
except ImportError:
    have_dbdreader = True
import glob
import logging
import numpy as np
import os
import re
import time
import xarray as xr
import xml.etree.ElementTree as ET
from collections.abc import Iterable

import pyglider.utils as utils


_log = logging.getLogger(__name__)


def binary_to_rawnc(indir, outdir, cacdir,
                    sensorlist, deploymentyaml,
                    incremental=True, scisuffix='EBD', glidersuffix='DBD'):
    """
    Convert slocum binary data (*.ebd/*.dbd) to raw netcdf files.

    Parameters
    ----------
    indir : str
        Directory with the raw ``*.ebd`` (science) and ``*.dbd`` (flight) files.
        These usually come from ``card_offload/Science/SENTLOGS`` or
        ``card_offload/Science/LOGS``, and ``card_offload/Main_board/SENTLOGS``
        and ``card_offload/Main_board/LOGS`. Recommend ``binary``

    outdir : str
        Directory to write the matching ``*.ebd.nc`` and ``*.dbd.nc`` files.
        Recommend ``rawnc``.

    cacdir : str
        Directory where the cached CAC sensor lists are kept.  These
        lists are often in directories like ``../Main_board/STATE/CACHE/``
        and ``../Science/STATE/CACHE/``, and the files in these directories
        should be copied to this directory by the user.  Recommend ``cac``

    sensorlist : str
        Text file with sensor list for this glider.  This filters the many
        sensors on the slocum gliders to just the ones listed here.  The file
        is text, comments deelineated as ``# a comment`` and a new entry
        on each line.

    deploymentyaml : str
        YAML text file with deployment information for this glider.

    incremental : bool, optional
        If *True* (default), only netcdf files that are older than the
        binary files are re-parsed.

    Returns
    -------
    status : bool
        *True* success.

    Notes
    -----
    This process can be slow for many files.
    """
    d = indir + '*.' + scisuffix
    filesScience = glob.glob(d)
    filesScience.sort()

    d = indir + '*.' + glidersuffix
    filesMain = glob.glob(d)
    filesMain.sort()

    keys = parse_filter_file(sensorlist)

    try:
        os.mkdir(outdir)
    except FileExistsError:
        pass

    todo = [filesMain, filesScience]
    sts = ['Flight', 'Science']
    for files, st in zip(todo, sts):
        _log.info(f'Working on {st}')
        # translate the flight files (sbd or dbd):
        deployment_ind_flight = 0
        badfiles = []
        for ind, filen in enumerate(files):
            # sometimes there is no science file for a flight file, so
            # we need to make sure the files match...
            try:
                fmeta, _ = dbd_get_meta(filen, cachedir=cacdir)
                path, ext = os.path.splitext(filen)
                fncname = (fmeta['the8x3_filename'] + '.' +
                           fmeta['filename_extension'] + '.nc')
                fullfncname = outdir + '/' + fncname

                ncfilesexist = os.path.isfile(fullfncname)
                if incremental and ncfilesexist:
                    ncfilesold = (os.path.getmtime(filen) >=
                                  os.path.getmtime(fullfncname))
                else:
                    ncfilesold = True
                if ncfilesold:
                    fdata, fmeta = dbd_to_dict(
                        filen, cacdir, keys=keys)
                    # need a unique index that increases monotonically, and is
                    # different for each file (we can't use time because it is
                    # not necessarily monotonic):
                    deployment_ind_flight = (int(fmeta['the8x3_filename']) * 1.0e6)
                    ds, deployment_ind_flight = datameta_to_nc(
                        fdata, fmeta, outdir=outdir, name=fncname,
                        deployment_ind=deployment_ind_flight)
            except Exception as e:
                badfiles += [files[ind]]
                _log.warning('Could not do parsing for %s', files[ind])
                _log.warning('%s', e)

        if len(badfiles) > 0:
            _log.warning('Some files could not be parsed:')
            for fn in badfiles:
                _log.warning('%s', fn)

    _log.info('All done!')


def _decode_sensor_info(dfh, meta):
    """
    Helper to decode the sensor list.

    dfh must be a filehandle because we want to be able to say where we stopped
    in file.
    """

    nsensors_total = int(meta['total_num_sensors'])
    nsensors_used = int(meta['sensors_per_cycle'])
    activeSensorList = [{} for i in range(nsensors_used)]
    outlines = []
    sensorInfo = {}
    for i in range(nsensors_total):
        line = dfh.readline().decode('utf-8')
        if line.split(':')[0] != 's':
            raise ValueError('Failed to parse sensor info')
        splitLine = [string.strip() for string in line.split(' ')[1:]
                     if string and not string.isspace()]
        sensorInfo[splitLine[-2]] = splitLine
        if splitLine[0] == 'T':
            activeSensorList[int(splitLine[2])] = {
                'name': splitLine[-2], 'unit': splitLine[-1],
                'bits': splitLine[-3]}
        outlines = outlines + [line]

    bindatafilepos = dfh.tell()  # keep this for seeking

    return activeSensorList, sensorInfo, outlines, bindatafilepos


def _get_cached_sensorlist(cachedir, meta):
    """
    Helper to get the sensor list from a file in the cache
    """
    fname0 = cachedir + '/' + meta['sensor_list_crc'].upper() + '.CAC'
    dd = glob.glob(cachedir + '/*')
    found = False
    for d in dd:
        if (os.path.split(d)[1].upper() ==
                os.path.split(fname0)[1].upper()):
            found = True
            break
    if not found:
        raise FileNotFoundError(f'Could not find {fname0}')

    with open(d, 'rb') as dfh:
        activeSensorList, sensorInfo, outlines, bindatafilepos = \
                _decode_sensor_info(dfh, meta)

    return activeSensorList, sensorInfo


def _make_cache(outlines, cachedir, meta):
    """
    Helper to make a cache file if one doesn't exist.
    """
    try:
        os.mkdir(cachedir)
    except FileExistsError:
        pass

    fname = cachedir + '/' + meta['sensor_list_crc'] + '.CAC'
    with open(fname, 'w') as dfh:
        for line in outlines:
            dfh.write(line)


def dbd_get_meta(filename, cachedir):
    """
    Get metadata from a dinkum binary file.

    Parameters
    ----------

    filename : str
        filename of the dinkum binary file (i.e. *.dbd, *.ebd)

    cachedir : str
        Directory where the cached CAC sensor lists are kept.  These
        lists are often in directories like ``../Main_board/STATE/CACHE/``.
        These should be copied somewhere locally.  Recommend ``./cac/``.

    Returns
    -------
    meta : dict
        Dictionary of the meta data for this dinkum binary file.

    """

    meta = {}

    with open(filename, 'rb') as dfh:
        meta['num_ascii_tags'] = 99  # read the first 99 lines
        while (len(meta) < int(meta['num_ascii_tags'])):
            line = dfh.readline().decode('utf-8')
            meta[line.split(':')[0]] = line.split(':')[1].strip()
        if len(meta) != int(meta['num_ascii_tags']):
            raise ValueError('Did not find expected number of tags')
        bindatafilepos = dfh.tell()
        localcache = False
        # if the sensor data is here, we need to read it, even though we
        # will use the cache
        if ('sensor_list_factored' in meta and
                not int(meta['sensor_list_factored'])):
            localcache = True
            activeSensorList, sensorInfo, outlines, bindatafilepos = \
                _decode_sensor_info(dfh, meta)

        # read the cache first.  If its not there, try to make one....
        try:
            activeSensorList, sensorInfo = \
                _get_cached_sensorlist(cachedir, meta)
        except FileNotFoundError:
            if localcache:
                _log.info('No cache file found; trying to create one')
                _make_cache(outlines, cachedir, meta)
            else:
                raise FileNotFoundError(
                    'No active sensor list found for crc ',
                    f'{meta["sensor_list_crc"]}. These are often found in ',
                    'offloaddir/Science/STATE/CACHE/ or ',
                    'offloaddir/Main_board/STATE/CACHE/. ',
                    f'Copy those locally into {cachedir}')
        meta['activeSensorList'] = activeSensorList
        # get the file's timestamp...
        meta['_dbdfiletimestamp'] = os.path.getmtime(filename)

    return meta, bindatafilepos


def dbd_to_dict(dinkum_file, cachedir, keys=None):
    """
    Translate a dinkum binary file to a dictionary of data and meta values.

    Parameters
    ----------
    dinkum_file : dbd file name (full path)
        These are the raw data from the glider, either offloaded from a card
        or from the dockserver.

    cachedir : str
        Directory where the cached CAC sensor lists are kept.  These
        lists are often in directories like ``../Main_board/STATE/CACHE/``.
        These should be copied somewhere locally.  Recommend ``./cac/``.

    keys : list of str
        list of sensor names to include in the *data* dictionary.  This
        allows us to make the dictionaries more compact and not have
        all the redundant sensor info.  If a single string then keys is a
        file name and passed to  `~.slocum.parse_filter_file` to get the list
        of keys.

    Returns
    -------
    data : dict
        dictionary of all the data with sensor names as keys, filtered
        according to the *keys* kwarg.

    meta : dict
        dictionary of all the meta data in the file.

    """
    # Parse ascii header - read in the metadata.
    data = []
    DINKUMCHUNKSIZE = int(3e4)  # how much data to pre-allocate

    if isinstance(keys, str):
        keys = parse_filter_file(keys)

    meta, bindatafilepos = dbd_get_meta(dinkum_file, cachedir)
    activeSensorList = meta['activeSensorList']
    dfh = open(dinkum_file, 'rb')
    # ------------------------------------------
    # All subsequent lines are in binary format.
    # Grab the seek pos and use that for a bookmark.
    # ------------------------------------------
    # offset for number of characters already read in.
    _log.debug('reading file from %d', bindatafilepos * 8)
    binaryData = bitstring.BitStream(dfh, offset=bindatafilepos * 8)
    # First there's the s,a,2byte int, 4 byte float, 8 byte double.
    # sometimes the endianess seems to get swapped.
    # ref_tuple = ['s', 'a', 4660, 123.456, 123456789.12345]
    diag_header = binaryData.readlist(['bits:8', 'bits:8'])
    diag_header[0] = chr(int(diag_header[0].hex, 16))
    diag_header[1] = chr(int(diag_header[1].hex, 16))
    if not (diag_header[0] == 's') and (diag_header[1] == 'a'):
        _log.warning("character failure: %s != 's', 'a'", diag_header)
        raise ValueError('Diagnostic header check failed.')

    endian = 'be'
    data = binaryData.read(f'uint{endian}:16')
    _log.debug('Checking endianness %s == 4660 or 13330', data)
    if data == 4660:
        pass
    elif data == 13330:
        endian = 'le'
    else:
        _log.warning("integer failure: %s != 4660", data)
        raise ValueError("Diagnostic header check failed.")
    _log.debug('Endianness is %s', endian)

    data = binaryData.read(f'float{endian}:32')
    if not np.allclose(data, 123.456):
        _log.warning("float32 failure: %s != 123.456", data)
        raise ValueError("Diagnostic header check failed.")

    data = binaryData.read(f'float{endian}:64')
    if not np.allclose(data, 123456789.12345):
        _log.warning("float64 failure: %s != 123456789.12345", data)
        raise ValueError("Diagnostic header check failed.")
    _log.debug('Diagnostic check passed.  Endian is %s', endian)

    nsensors = int(meta['sensors_per_cycle'])
    currentValues = np.zeros(int(meta['sensors_per_cycle'])) + np.nan
    data = np.zeros((DINKUMCHUNKSIZE, nsensors)) + np.nan
    # Then there's a data cycle with every sensor marked as updated, giving
    # us our initial values.
    # 01 means updated with 'same value', 10 means updated with a new value,
    # 11 is reserved, 00 is not updated.
    # This character leads off each byte cycle.
    frameCheck = binaryData.read('bytes:1').decode("utf-8")
    updatedCode = ['00'] * int(meta['sensors_per_cycle'])

    # Data cycle begins now.
    # Cycle tag is a ascii 'd' character. Then
    # state_bytes_per_cycle * state_bytes (2bits per sensor) of state bytes.
    # Then data for each updated sensor as per the state bytes.
    # Then zeroes until the last byte is completed, should they be necessary.
    _log.info('Parsing binary data')
    proctimestart = time.time()
    ndata = 0
    while frameCheck == 'd':
        for i in range(int(meta['sensors_per_cycle'])):
            updatedCode[i] = binaryData.read('bin:2')
        # burn off any remaining bits to get to the first full bit.
        binaryData.bytealign()
        for i, code in enumerate(updatedCode):
            if code == '00':  # No new value
                currentValues[i] = np.nan
            elif code == '01':  # Same value as before.
                continue
            elif code == '10':  # New value.
                if int(activeSensorList[i]['bits']) in [4, 8]:
                    currentValues[i] = binaryData.read(
                        f'float{endian}:' +
                        str(int(activeSensorList[i]['bits']) * 8))
                elif int(activeSensorList[i]['bits']) in [1, 2]:
                    currentValues[i] = binaryData.read(
                        f'uint{endian}:' +
                        str(int(activeSensorList[i]['bits']) * 8))
                else:
                    raise ValueError('Bad bits')
            else:
                raise ValueError(('Unrecognizable code in data cycle. ',
                                  'Parsing failed'))
        data[ndata] = currentValues
        binaryData.bytealign()

        # We've arrived at the next line.
        try:
            d = binaryData.peek('bytes:1').decode('utf-8')
        except bitstring.ReadError:
            _log.debug('position at end of stream %d',
                       binaryData.pos + 8 * bindatafilepos)
            _log.warning('End of file reached without termination char')
            d = 'X'
        if d == 'd':
            frameCheck = binaryData.read('bytes:1').decode('utf-8')
            ndata += 1
            if ndata % DINKUMCHUNKSIZE == 0:
                # need to allocate more data!
                data = np.concatenate(
                    (data, np.nan + np.zeros((DINKUMCHUNKSIZE, nsensors))),
                    axis=0)
        elif d == 'X':
            # End of file cycle tag. We made it through.
            # throw out pre-allocated data we didn't use...
            data = data[:ndata]
            break
        else:
            raise ValueError(f'Parsing failed at {binaryData.bytepos}. ',
                             f'Got {d} expected d or X')

    proctimeend = time.time()
    _log.info(('%s lines of data read from %s, data rate of %s rows '
               'per second') % (len(data), dinkum_file,
                                len(data) / (proctimeend - proctimestart)))
    dfh.close()

    _log.info('Putting data into dictionary')
    ddict = dict()

    # deal 2-D array into a dictionary...  Only keep keys we want...
    for n, key in enumerate(meta['activeSensorList']):
        if keys is None or key['name'] in keys:
            ddict[key['name']] = data[:, n]

    return ddict, meta


def add_times_flight_sci(fdata, sdata=None):
    """
    Add the time from the flight data to the science data structure.

    Parameters
    ----------
    fdata, sdata : dict
        data dictionaries from ``dbd_to_dict``.  If sdata = None assume
        no sdata file given.

    Returns
    -------
    fdata, sdata : dict
        as input, but with 'm_present_time_sci' added to sdata and
        'sci_m_present_time_fixed' added to fdata.
    """

    # Basically throw out the leading flight timestamp if its lag threshhold
    # is too high.

    uniqueTimes, uniqueTimeIndices = np.unique(fdata['m_present_time'],
                                               return_index=True)
    if len(uniqueTimes) != len(fdata['m_present_time']):
        # Correct the duplicates in the flight timestamps..
        _log.warning('Duplicate flight entries detected.')
    # Fix common problems with science data set.
    fdata['sci_m_present_time_fixed'] = (fdata['sci_m_present_time'] +
                                         fdata['m_science_clothesline_lag'])
    # There are some nans in the sci_m_present_time_fixed set.
    # We need to interpolate them.

    # Interpolate the nans out of sci_m_present_time.
    good = ~np.isnan(fdata['sci_m_present_time_fixed'])
    bad = ~good
    if not np.all(bad):
        fdata['sci_m_present_time_fixed'][bad] = np.interp(
            fdata['m_present_time'][bad], fdata['m_present_time'][good],
            fdata['sci_m_present_time_fixed'][good])

    if sdata is not None:
        lag_threshhold = np.nanmax(fdata['u_max_clothesline_lag_for_consci'])
        _log.info('lag_threshhold %f', lag_threshhold)
        # Number of seconds the computers can be apart before we stop believing
        # them. Crucial for times the science computer is stopped.
        # Calculate the equivalent flight computer times for each science
        # timestamp, given the common times we know about and have verified as
        # unique.
        # If you are inclined to believe the flight timestamps. Throw out and
        # interpolate if the time is greater than the set lag threshhold:

        # With the recommended fix from Kerfoot's lab and D. Pingal
        # Some of the m_present_times are nan-ed...
        good = np.logical_and(
            (fdata['m_science_clothesline_lag'][uniqueTimeIndices] <
                lag_threshhold),
            np.isfinite(fdata['m_science_clothesline_lag'][uniqueTimeIndices]))

        if not np.all(~good):
            tf = fdata['sci_m_present_time_fixed'][uniqueTimeIndices[good]]
            pt = fdata['m_present_time'][uniqueTimeIndices[good]]

            sdata['m_present_time_sci'] = np.interp(
                sdata['sci_m_present_time'], tf, pt, np.nan, np.nan)
        else:
            sdata['m_present_time_sci'] = np.nan * sdata['sci_m_present_time']

    return fdata, sdata


def parse_filter_file(filter_file):
    keys = []
    with open(filter_file) as fin:
        for li in fin:
            if not li[0] in ['#']:
                lis = li.split(' ')
                if len(lis) == 1:
                    key = lis[0].rstrip()
                    if len(key) > 0:
                        keys += [key]
    return keys


def datameta_to_nc(data, meta, outdir=None, name=None, check_exists=False,
                   deployment_ind=0):
    """
    Convert a raw dinkum data and meta dict to a netcdf file.

    Parameters
    ----------
    data, meta : dict
        data, meta are a pair of dicts returned from `dbd_to_dict`.

    outdir : str
        directory where the netcdf file will be written.

    name : str or None
        name of the file (including extension).  If None the *name* of
        the file will be taken from ``meta['full_filename']``.

    check_exists : bool
        If the netcdf file exists, and meta['_dbdfiletimestamp'] is older than
        the netcdf file modified time, don't remake the netcdf file.  Default
        is False, and the netcdf file will be overwritten.
    """

    if name is None:
        name = meta['full_filename'] + '.' + meta['filename_extension'] + '.nc'
    if outdir is None:
        outdir = './'
    outname = outdir + '/' + name

    if check_exists:
        if (os.path.isfile(outname) and
                (os.path.getmtime(outname) > meta['_dbdfiletimestamp'])):
            _log.info('%s already exists and is newer than timestamp in meta',
                      outname)
            return None

    ds = xr.Dataset()
    if 'm_present_time' in data.keys():
        time = data['m_present_time']
    if 'sci_m_present_time' in data.keys():
        time = data['sci_m_present_time']

    index = np.arange(len(time)) + deployment_ind
    # this gets passed to the next file...
    ds['_ind'] = (('_ind'), index)

    ds['time'] = (('_ind'), time)
    for key in data.keys():
        ds[key] = (('_ind'), data[key])
        # try and find the unit for this....
        for sensor in meta['activeSensorList']:
            if sensor['name'] == key:
                ds[key].attrs['unit'] = sensor['unit']
                break
    for key in meta.keys():
        if key != 'activeSensorList':
            ds.attrs[key] = meta[key]
        # make a long string for activeSensorList:
        listst = ''
        for sensor in meta['activeSensorList']:
            listst += '%s' % sensor
            listst += '\n'
        ds.attrs['activeSensorList'] = listst

    ds.attrs['_processing'] = __name__ + ' python library'
    ds.attrs['Conventions'] = 'None'

    # trim data that has time==0
    ind = np.where(ds.time > 1e4)[0]
    # _log.debug(f'{ds}, {ds.time}, {len(ind)}')
    ds = ds.isel(_ind=ind)
    ds['_ind'] = np.arange(len(ds.time)) + deployment_ind
    if len(ds['_ind'].values) > 1:
        deployment_ind = ds['_ind'].values[-1]

    _log.info(f'Writing! {deployment_ind}')
    ds.to_netcdf(outname, 'w')
    _log.info(f'Wrote:, {outname}')
    return ds, deployment_ind


def merge_rawnc(indir, outdir, deploymentyaml,
                scisuffix='EBD', glidersuffix='DBD'):
    """
    Merge all the raw netcdf files in indir.  These are meant to be
    the raw flight and science files from the slocum.

    Parameters
    ----------
    indir : str
        Directory where the raw ``*.ebd.nc`` and ``*.dbd.nc`` files are.
        Recommend: ``./rawnc``

    outdir : str
        Directory where merged raw netcdf files will be put. Recommend:
        ``./rawnc/``.  Note that the netcdf files will be named following
        the data in *deploymentyaml*:
        ``glider_nameglider_serial-YYYYmmddTHHMM-rawebd.nc`` and
        ``...rawdbd.nc``.

    deploymentyaml : str
        YAML text file with deployment information for this glider.

    """

    scisuffix = scisuffix.lower()
    glidersuffix = glidersuffix.lower()
    deployment = utils._get_deployment(deploymentyaml)

    metadata = deployment['metadata']
    id = metadata['glider_name'] + metadata['glider_serial']

    # different missions get a different number and they may not merge
    # smoothly, hence the different files made here.

    # first weed out singleton files.  These cause merge problems...
    d = indir + '/*.' + scisuffix + '.nc'
    fin = glob.glob(d)
    for f in fin:
        with xr.open_dataset(f) as ds:
            if len(ds._ind) < 2:
                bad = True
            else:
                bad = False
        if bad:
            os.rename(f, f+'.singleton')
            fglider = f
            try:
                fglider = fglider.replace(scisuffix, glidersuffix)
                os.rename(fglider, fglider+'.singleton')
            except FileNotFoundError:
                pass

    fin = glob.glob(indir + '/*.' + glidersuffix + '.nc')
    with xr.open_mfdataset(fin, decode_times=False, lock=False) as ds:
        outnebd = outdir + '/' + id + 'rawdbd.nc'
        ds = ds.sortby('time')
        ds['_ind'] = np.arange(len(ds.time))
        ds.to_netcdf(outnebd, 'w')

    fin = glob.glob(indir + '/*.' + scisuffix + '.nc')

    with xr.open_mfdataset(fin, decode_times=False, lock=False) as ds:
        outnebd = outdir + '/' + id + 'rawebd.nc'
        ds = ds.sortby('time')
        ds['_ind'] = np.arange(len(ds.time))
        ds.to_netcdf(outnebd, 'w')


def raw_to_timeseries(indir, outdir, deploymentyaml, *,
                      profile_filt_time=100, profile_min_time=300):
    """
    Parameters
    ----------
    indir : string
        Directory with raw netcdf files.
    outdir : string
        Directory to put the merged timeseries files.
    profile_filt_time : float
        time in seconds over which to smooth the pressure time series for
        finding up and down profiles (note, doesn't filter the data that is
        saved)
    profile_min_time : float
        minimum time to consider a profile an actual profile (seconds)
    Returns
    -------
    outname : string
        name of the new merged netcdf file.
    """

    deployment = utils._get_deployment(deploymentyaml)
    metadata = deployment['metadata']
    ncvar = deployment['netcdf_variables']
    device_data = deployment['glider_devices']
    thenames = list(ncvar.keys())
    thenames.remove('time')

    id = metadata['glider_name'] + metadata['glider_serial']

    id0 = None
    ebdn = indir + '/' + id + 'rawebd.nc'
    dbdn = indir + '/' + id + 'rawdbd.nc'
    _log.debug(f'{ebdn}, {dbdn}')
    if not os.path.exists(ebdn) or not os.path.exists(dbdn):
        raise FileNotFoundError('Could not find %s and %s', ebdn, dbdn)

    _log.info(f'Opening: {ebdn}, {dbdn}')
    ebd = xr.open_dataset(ebdn, decode_times=False)
    dbd = xr.open_dataset(dbdn, decode_times=False)
    _log.debug(f'DBD, {dbd}, {dbd.m_depth}')
    if len(ebd.time) <= 2:
        raise RuntimeError('Science file has no data')
    # build a new data set based on info in `deployment.`
    # We will use ebd.m_present_time as the interpolant if the
    # variable is in dbd.
    ds = xr.Dataset()
    attr = {}
    name = 'time'
    for atts in ncvar[name].keys():
        if atts != 'coordinates':
            attr[atts] = ncvar[name][atts]
    ds[name] = (('time'), ebd[name].values, attr)
    for name in thenames:
        _log.info('working on %s', name)
        if 'method' in ncvar[name].keys():
            continue
        # variables that are in the data set or can be interpolated from it
        if 'conversion' in ncvar[name].keys():
            convert = getattr(utils, ncvar[name]['conversion'])
        else:
            convert = utils._passthrough
        sensorname = ncvar[name]['source']
        _log.info('names: %s %s', name, sensorname)
        if sensorname in ebd.keys():
            _log.debug('EBD sensorname %s', sensorname)
            val = ebd[sensorname]
            val = utils._zero_screen(val)
    #        val[val==0] = np.nan
            val = convert(val)
        else:
            _log.debug('DBD sensorname %s', sensorname)
            val = convert(dbd[sensorname])
            val = _dbd2ebd(dbd, ds, val)
            ncvar['method'] = 'linear fill'
        # Check for NMEA to degrees conversion
        if 'source_units' in ncvar[name] and 'units' in ncvar[name]:
            if ncvar[name]['source_units'] == 'NMEA':
                if 'degrees' in ncvar[name]['units']:
                    _log.debug(f'Converting {name} NMEA to degrees')
                    val.data = utils.nmea2deg(val.data)
        # make the attributes:
        ncvar[name].pop('coordinates', None)
        attrs = ncvar[name]
        attrs = utils.fill_required_attrs(attrs)
        ds[name] = (('time'), val.data, attrs)

    _log.debug(f'HERE, {ds}')
    _log.debug(f'HERE, {ds.pressure[0:100]}')
    # some derived variables:
    # trim bad times...

    ds = utils.get_glider_depth(ds)
    ds = utils.get_distance_over_ground(ds)

    ds = utils.get_derived_eos_raw(ds)
    ds = ds.assign_coords(longitude=ds.longitude)
    ds = ds.assign_coords(latitude=ds.latitude)
    ds = ds.assign_coords(depth=ds.depth)

    ds['time'] = (ds.time.values.astype('timedelta64[s]') +
                  np.datetime64('1970-01-01T00:00:00')).astype('datetime64[ns]')
    _log.info(f'utils.fill_metadata: {device_data}')
    ds = utils.fill_metadata(ds, deployment['metadata'], device_data)
    start = ds['time'].values[0]
    end = ds['time'].values[-1]

    ds.attrs['deployment_start'] = str(start)
    ds.attrs['deployment_end'] = str(end)
    _log.debug(ds.depth.values[:100])
    _log.debug(ds.depth.values[2000:2100])
    ds = utils.get_profiles_new(
        ds, filt_time=profile_filt_time, profile_min_time=profile_min_time)
    _log.debug(ds.depth.values[:100])
    _log.debug(ds.depth.values[2000:2100])

    try:
        os.mkdir(outdir)
    except:
        pass
    outname = (outdir + '/' + ds.attrs['deployment_name'] + '.nc')
    _log.info('writing %s', outname)
    ds.to_netcdf(outname, 'w',
                 encoding={'time': {'units': 'seconds since 1970-01-01T00:00:00Z'}})
    if id0 is None:
        id0 = ds.attrs['deployment_name']

    return outname


def binary_to_timeseries(indir, cachedir, outdir, deploymentyaml, *,
                         search='*.[D|E]BD', fnamesuffix='',
                         time_base='sci_water_temp', profile_filt_time=100,
                         profile_min_time=300, maxgap=300,
                         replace_attrs=None):
    """
    Convert directly from binary files to netcdf timeseries file.  Requires
    dbdreader to be installed.

    Parameters
    ----------
    indir : string
        Directory with binary files from the glider.

    cachedir : string
        Directory with glider cache files (cac files)

    outdir : string
        Directory to put the merged timeseries files.

    deploymentyaml : str or list
        Name of YAML text file with deployment information for this glider.

        If a list, then the YAML files are read in order, and any top-level dictionaries
        are overwritten from the previous YAMLs.  The advantage of this is that it allows
        metadata that is common to multiple ways of processing the data come from the
        first file, and then subsequent files change "netcdf_variables" if desired.

    profile_filt_time : float or None
        time in seconds over which to smooth the pressure time series for
        finding up and down profiles (note, doesn't filter the data that is
        saved).  If None, then do not find profiles.

    profile_min_time : float or None
        minimum time to consider a profile an actual profile (seconds).  If None,
        then do not find profiles.

    maxgap : float
        Longest gap in seconds to interpolate over when matching instrument
        timeseries.

    replace_attrs : dict or None
        replace global attributes in the metadata after reading the metadata
        file in.  Helpful when processing runs with only a couple things that
        change.

    Returns
    -------
    outname : string
        name of the new merged netcdf file.
    """

    if not have_dbdreader:
        raise ImportError('Cannot import dbdreader')

    deployment = utils._get_deployment(deploymentyaml)
    if replace_attrs:
        for att in replace_attrs:
            deployment['metadata'][att] = replace_attrs[att]

    ncvar = deployment['netcdf_variables']
    device_data = deployment['glider_devices']
    thenames = list(ncvar.keys())
    thenames.remove('time')

    # get the dbd file
    _log.info(f'{indir}/{search}')
    dbd = dbdreader.MultiDBD(pattern=f'{indir}/{search}',
                             cacheDir=cachedir)

    # build a new data set based on info in `deployment.`
    # We will use ebd.m_present_time as the interpolant if the
    # variable is in dbd.
    ds = xr.Dataset()
    attr = {}
    name = 'time'
    for atts in ncvar[name].keys():
        if (atts != 'coordinates') & (atts != 'units') & (atts != 'calendar'):
            attr[atts] = ncvar[name][atts]
    sensors = [time_base]

    for nn, name in enumerate(thenames):
        sensorname = ncvar[name]['source']
        if not sensorname == time_base:
            sensors.append(sensorname)
        else:
            baseind = nn

    # get the data, with `time_base` as the time source that
    # all other variables are synced to:
    data = list(dbd.get_sync(*sensors))
    # get the time:
    time = data.pop(0)
    ds['time'] = (('time'), time, attr)
    ds['latitude'] = (('time'), np.zeros(len(time)))
    ds['longitude'] = (('time'), np.zeros(len(time)))
    # get the time_base data:
    basedata = data.pop(0)
    # slot the time_base variable into the right place in the
    # data list:
    data.insert(baseind, basedata)

    for nn, name in enumerate(thenames):
        _log.info('working on %s', name)
        if 'method' in ncvar[name].keys():
            continue
        # variables that are in the data set or can be interpolated from it
        if 'conversion' in ncvar[name].keys():
            convert = getattr(utils, ncvar[name]['conversion'])
        else:
            convert = utils._passthrough

        sensorname = ncvar[name]['source']
        _log.info('names: %s %s', name, sensorname)
        if sensorname in dbd.parameterNames['sci']:
            _log.debug('Sci sensorname %s', sensorname)
            val = data[nn]

            # interpolate only over those gaps that are smaller than 'maxgap'
            # in seconds
            _t, _ = dbd.get(ncvar[name]['source'])
            tg_ind = utils.find_gaps(_t, time, maxgap)
            val[tg_ind] = np.nan

            val = utils._zero_screen(val)
            val = convert(val)
        elif sensorname in dbd.parameterNames['eng']:
            _log.debug('Eng sensorname %s', sensorname)
            val = data[nn]
            val = convert(val)
            ncvar['method'] = 'linear fill'
        else:
            ValueError(f'{sensorname} not in science or eng parameter names')

        # make the attributes:
        ncvar[name]['coordinates'] = 'time'
        attrs = ncvar[name]
        attrs = utils.fill_required_attrs(attrs)
        ds[name] = (('time'), val, attrs)


    _log.info(f'Getting glider depths, {ds}')
    _log.debug(f'HERE, {ds.pressure[0:100]}')

    ds = utils.get_glider_depth(ds)
    ds = utils.get_distance_over_ground(ds)

    ds = utils.get_derived_eos_raw(ds)

    # screen out-of-range times; these won't convert:
    ds['time'] = ds.time.where((ds.time>0) & (ds.time<6.4e9), np.nan)
    # convert time to datetime64:
    ds['time'] = (ds.time*1e9).astype('datetime64[ns]')
    ds['time'].attrs = attr

    ds = utils.fill_metadata(ds, deployment['metadata'], device_data)

    start = ds.time.values[0]
    end = ds.time.values[0]
    _log.debug('Long')
    _log.debug(ds.longitude.values[-2000:])
    # make sure this is ISO readable....
    ds.attrs['deployment_start'] = str(start)[:18]
    ds.attrs['deployment_end'] = str(end)[:18]
    _log.debug(ds.depth.values[:100])
    _log.debug(ds.depth.values[2000:2100])

    if (profile_filt_time is not None) and (profile_min_time is not None):
        ds = utils.get_profiles_new(
            ds, filt_time=profile_filt_time, profile_min_time=profile_min_time)
    _log.debug(ds.depth.values[:100])
    _log.debug(ds.depth.values[2000:2100])

    try:
        os.mkdir(outdir)
    except:
        pass
    outname = (outdir + '/' + ds.attrs['deployment_name'] + fnamesuffix + '.nc')
    _log.info('writing %s', outname)
    # convert time back to float64 seconds for ERDDAP etc happiness, as they won't take ns
    # as a unit:
    ds.to_netcdf(outname, 'w',
                 encoding={'time': {'units': 'seconds since 1970-01-01T00:00:00Z',
                                    '_FillValue': np.nan,
                                    'dtype': 'float64'}})

    return outname


# alias:
raw_to_L1timeseries = raw_to_L0timeseries = raw_to_timeseries


def timeseries_get_profiles(inname, profile_filt_time=100,
                            profile_min_time=400):
    """
    Parameters
    ----------
    profile_filt_time : float
        how long a filter to apply to the pressure data in seconds

    profile_min_time : float
        how long a profile must last to be considered a proper profile (seconds)
    """
    with xr.open_dataset(inname) as ds:
        ds = utils.get_profiles_new(
            ds, filt_time=profile_filt_time, profile_min_time=profile_min_time)
    ds.to_netcdf(inname, mode='a')
    return inname


def _dbd2ebd(dbd, ds, val):
    """
    Helper to interpolate from dbd to ebd data stream
    """
    good = np.where(np.isfinite(val))[0]
    vout = ds.time * 0.0
    goodt = np.where(np.isfinite(ds.time))[0]
    if (len(goodt) > 1) and (len(good) > 1):
        _log.debug(f'GOOD, {goodt}, {good}')
        vout[goodt] = np.interp(
            ds.time[goodt].values, dbd.m_present_time.values[good], val[good].values)
    return vout


def parse_gliderState(fname):
    """
    Parse time, lat, and lon from a gliderstate file

    Parameters
    ----------
    fname : string or Path
        Location of the gliderState.xml file for the glider

    Returns
    -------
    dat : xarray
        xarray with fields time, lon, lat
    """

    data = {'lon': ('report', np.zeros(10000)),
            'lat': ('report', np.zeros(10000)),
            'time': ('report', np.zeros(10000,
                                        dtype='datetime64[s]'))}

    dat = xr.Dataset(data_vars=data)

    tree = ET.parse(fname)
    root = tree.getroot()
    nevents = 0
    for event in root.iter('report'):
        e = event.find('locations')
        for el in e.iter('valid_location'):
            lat = el.find('lat').text
            lon = el.find('lon').text
            time = el.find('time').text
            if time != 'unavailable':
                dat.time[nevents] = np.datetime64(time[:-5])
                dat.lon[nevents] = utils.nmea2deg(float(lon))
                dat.lat[nevents] = utils.nmea2deg(float(lat))
                nevents += 1
    dat = dat.isel(report=slice(nevents))
    return dat


def parse_logfiles(files):
    """
    Parse time, lat, lon, and amph_total from glider logfiles.

    Parameters
    ----------
    files : list of strings or Paths
        List of logfiles to parse.  Should be sorted.

    Returns
    -------
    out : xarray
        xarray data set with fields time, lon, lat, ampH indexed by surfacing.
        More could be added.
    """

    times = [''] * 10000
    gps = [''] * 10000
    amph = [''] * 10000
    relcharge = np.zeros(10000) * np.nan
    volts = np.zeros(10000) * np.nan

    ntimes = 0
    for fn in files:
        found_time = False

        with open(fn, 'r') as fin:
            for ll in fin:
                if 'Curr Time:' in ll:
                    times[ntimes] = ll
                    ntimes += 1
                    found_time = True
                if found_time and 'GPS Location' in ll:
                    gps[ntimes - 1] = ll
                if found_time and "sensor:m_coulomb_amphr_total" in ll:
                    amph[ntimes-1] = ll
                if ll.startswith('   sensor:m_lithium_battery_relative_charge'):
                        pattern = r'=(\d+\.\d+)'
                        match = re.search(pattern, ll)
                        relcharge[ntimes-1] = float(match.group(1))
                if ll.startswith('   sensor:m_battery(volts)='):
                        pattern = r'=(\d+\.\d+)'
                        match = re.search(pattern, ll)
                        volts[ntimes-1] = float(match.group(1))
    amph = amph[:ntimes]
    gps = gps[:ntimes]
    times = times[:ntimes]
    volts = volts[:ntimes]
    relcharge = relcharge[:ntimes]
    # now parse them
    out = xr.Dataset(
        coords={'time': ('surfacing', np.zeros(ntimes, dtype='datetime64[ns]'))})
    out['ampH'] = ('surfacing', np.zeros(ntimes) * np.nan)
    out['lon'] = ('surfacing', np.zeros(ntimes) * np.nan)
    out['lat'] = ('surfacing', np.zeros(ntimes) * np.nan)
    # these don't need to be parsed:
    out['volts'] = ('surfacing', volts[:ntimes])
    out['relcharge'] = ('surfacing', relcharge[:ntimes])

    for i in range(ntimes):
        timestring = times[i][11:-13]
        try:
            out['time'][i] = np.datetime64(
                datetime.strptime(timestring, '%a %b %d %H:%M:%S %Y'), 'ns')
            st = amph[i].index('=')
            en = amph[i][st:].index(' ') + st
            out['ampH'][i] = float(amph[i][(st+1):en])
            sp = gps[i].split(' ')
            out['lat'][i] = utils.nmea2deg(float(sp[3]))
            out['lon'][i] = utils.nmea2deg(float(sp[5]))
        except:
            pass

    return out


def parse_logfiles_maybe(files):
    """
    Parse time, lat, lon, and amph_total from glider logfiles.

    Parameters
    ----------
    files : list of strings or Paths
        List of logfiles to parse.  Should be sorted.

    Returns
    -------
    out : xarray
        xarray data set with fields time, lon, lat, ampH indexed by surfacing.
        More could be added.
    """

    times = [''] * 10000
    gps = [''] * 10000
    amph = [''] * 10000
    surfacereason = ''
    missionnum = [''] * 10000
    abortsegment = 0
    abortcause = 0

    ntimes = 0
    for fn in files:
        found_time = False

        with open(fn, 'r') as fin:
            for l in fin:
                if 'Curr Time:' in l:
                    times[ntimes] = l
                    ntimes += 1
                    found_time=True
                elif found_time and 'GPS Location' in l:
                    gps[ntimes - 1] = l
                elif found_time and "sensor:m_coulomb_amphr_total" in l:
                    amph[ntimes-1] = l
                elif found_time and "Because:" in l:
                    surfacereason = l
                elif found_time and "MissionNum" in l:
                    missionnum[ntimes-1] = l
                elif found_time and "abort segment:" in l:
                    abortsegment = l
                elif found_time and "abort cause:" in l:
                    abortcause = l

    amph = amph[:ntimes]
    gps = gps[:ntimes]
    times = times[:ntimes]
    missionnum = missionnum[:ntimes]

    # now parse them
    out = xr.Dataset(coords={'time': ('surfacing', np.zeros(ntimes, dtype='datetime64[ns]'))})
    out['ampH'] = ('surfacing', np.zeros(ntimes) * np.nan)
    out['lon'] = ('surfacing', np.zeros(ntimes) * np.nan)
    out['lat'] = ('surfacing', np.zeros(ntimes) * np.nan)
    out['missionnum'] = ('surfacing', np.zeros(ntimes) * np.nan)
    out.attrs['surfacereason'] = surfacereason
    # ABORT HISTORY: last abort segment: hal_1002-2024-183-0-0 (0171.0000)
    out.attrs['abortsegment'] = float(abortsegment[-11:-2])
    out.attrs['abortcause'] = abortcause

    for i in range(ntimes):
        timestring = times[i][11:-13]
        out['time'][i] = np.datetime64(datetime.strptime(timestring, '%a %b %d %H:%M:%S %Y'), 'ns')
        try:
            if '=' in amph[i]:
                st = amph[i].index('=')
                en = amph[i][st:].index(' ') + st
                out['ampH'][i] = float(amph[i][(st+1):en])

            #        GPS Location:  4912.737 N -12357.253 E measured    110.757 secs ago
            sp = gps[i].split()
            out['lat'][i] = utils.nmea2deg(float(sp[2]))
            out['lon'][i] = utils.nmea2deg(float(sp[4]))

            # MissionName:calvert.mi MissionNum:hal_1002-2024-183-4-41 (0175.0041)
            if len(missionnum[i]) > 12:
                out['missionnum'][i] = float(missionnum[i][-11:-2])
        except:
            pass

    return out





__all__ = ['binary_to_rawnc', 'merge_rawnc', 'raw_to_timeseries',
           'parse_gliderState', 'parse_logfiles']
