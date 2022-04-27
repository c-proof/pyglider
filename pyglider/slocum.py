# -*- coding: utf-8 -*-
from pyglider import bitstring
import datetime
import glob
import itertools
import logging
from math import floor, fmod
import numpy as np
import os
import re
import time
import xarray as xr
import yaml
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
        Directory with the raw *.ebd (science) and *.dbd (flight) files.
        These usually come from ``card_offload/Science/SENTLOGS`` or
        ``card_offload/Science/LOGS``, and ``card_offload/Main_board/SENTLOGS`
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

    if len(filesMain) < 0 or len(filesScience) < 0:
        raise FileNotFoundError(
            'No science or flight files found in %s' % indir)

    keys = parse_filter_file(sensorlist)

    try:
        os.mkdir(outdir)
    except FileExistsError:
        pass

    deployment_ind_sci = 0
    deployment_ind_flight = 0
    badfiles = []
    for ind in range(len(filesMain)):
        # sometimes there is no science file for a flight file, so
        # we need to make sure the files match...
        try:
            fmeta, _ = dbd_get_meta(filesMain[ind], cachedir=cacdir)
            path, ext =  os.path.splitext(filesMain[ind])
            sciname = filesMain[ind][:-3] + scisuffix
            fncname = (fmeta['the8x3_filename'] + '.' +
                       fmeta['filename_extension'] + '.nc')
            fullfncname = outdir + '/' + fncname

            if os.path.isfile(sciname):
                smeta, _ = dbd_get_meta(sciname, cachedir=cacdir)
                sncname = (smeta['the8x3_filename'] + '.' +
                           smeta['filename_extension'] + '.nc')
                fullsncname = outdir + '/' + sncname

                _log.info('Working on  %s %s',
                    os.path.basename(filesMain[ind]),
                    os.path.basename(sciname))

                ncfilesexist = (os.path.isfile(fullsncname) and
                                os.path.isfile(fullfncname))
                if incremental and ncfilesexist:
                    ncfilesold = ((os.path.getmtime(sciname) >=
                                   os.path.getmtime(fullsncname)) and
                                  (os.path.getmtime(filesMain[ind]) >=
                                   os.path.getmtime(fullfncname)))
                else:
                    ncfilesold = True
                if ncfilesold:
                    # save these in case they get corrupted below...
                    dis = deployment_ind_sci
                    dif = deployment_ind_flight
                    _log.info(f'sciind: {dis}; gliderind: {dif}')
                    if 1:
                        sdata, smeta = dbd_to_dict(sciname, cacdir, keys=keys)
                        fdata, fmeta = dbd_to_dict(filesMain[ind], cacdir,
                                keys=keys)
                        deployment_ind_sci = int(smeta['the8x3_filename'])*1.0e4
                        ds, deployment_ind_sci = datameta_to_nc(sdata,
                            smeta, outdir=outdir,
                            name=sncname, deployment_ind=deployment_ind_sci)
                        deployment_ind_flight = (int(smeta['the8x3_filename']) *
                                                 1.0e4)
                        ds, deployment_ind_flight = datameta_to_nc(fdata, fmeta,
                            outdir=outdir, name=fncname,
                            deployment_ind=deployment_ind_flight)
                    else:
                        deployment_ind_sci = dis
                        deployment_ind_flight = dif
                        _log.warning('Could not decode %s', filesScience[ind])
                else:
                    _log.info('skipping %s', sciname)
            else:
                _log.info('No science file found for %s', filesMain[ind])

        except:
            badfiles += [filesMain[ind]]
            _log.warning('Could not do parsing for %s', filesMain[ind])
        _log.info('')

    if len(badfiles) > 0:
        _log.warning('Some files could not be parsed:')
        for fn in badfiles:
            _log.warning('%s', fn)

    _log.info('All done!')



def _check_diag_header(diag_tuple):
    # diagnostic values should be
    # ['s', 'a', 4660, 123.45600128173828, 123456789.12345] # 4660 is 0x1234
    ref_tuple = ['s', 'a', 4660, 123.456, 123456789.12345]
    for i in range(3):
        if ref_tuple[i] != diag_tuple[i]:
            _log.warning('character or int failure: %s', diag_tuple)
            return False

    if ((abs(ref_tuple[3] - diag_tuple[3]) > .0001) or
        (abs(ref_tuple[4] - diag_tuple[4]) > .0001)):
        _log.warning('floating point failure')
        return False
    return True


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
                raise FileNotFoundError(('No active sensor list found for crc ',
                    '{}. These are often found in ',
                    'offloaddir/Science/STATE/CACHE/ or ',
                    'offloaddir/Main_board/STATE/CACHE/. ',
                    'Copy those locally into {}'
                    ).format(meta['sensor_list_crc'], cachedir))
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
        file name and passed to  `dinkum.parse_filter_file` to get the list
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
    binaryData = bitstring.BitStream(dfh, offset=bindatafilepos * 8)
    # First there's the s,a,2byte int, 4 byte float, 8 byte double.
    diag_header = binaryData.readlist(['bits:8', 'bits:8',
                                       'uint:16', 'float:32', 'float:64'])
    diag_header[0] = chr(int(diag_header[0].hex, 16))
    diag_header[1] = chr(int(diag_header[1].hex, 16))
    if not _check_diag_header(diag_header):
        raise ValueError('Diagnostic header check failed.')

    nsensors = int(meta['sensors_per_cycle'])
    currentValues = np.zeros(int(meta['sensors_per_cycle'])) + np.NaN
    data = np.zeros((DINKUMCHUNKSIZE, nsensors)) + np.NaN
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
                currentValues[i] = np.NaN
            elif code == '01':  # Same value as before.
                continue
            elif code == '10':  # New value.
                if int(activeSensorList[i]['bits']) in [4, 8]:
                    currentValues[i] = binaryData.read(
                        'float:' + str(int(activeSensorList[i]['bits']) * 8))
                    # print currentValues[i], activeSensorList[i]['name']
                elif int(activeSensorList[i]['bits']) in [1, 2]:
                    currentValues[i] = binaryData.read(
                        'uint:' + str(int(activeSensorList[i]['bits']) * 8))
                else:
                    raise ValueError('Bad bits')
            else:
                raise ValueError(('Unrecognizable code in data cycle. ',
                                  'Parsing failed'))
        data[ndata] = currentValues
        binaryData.bytealign()
        # We've arrived at the next line.
        d = binaryData.peek('bytes:1').decode('utf-8')
        if d == 'd':
            frameCheck = binaryData.read('bytes:1').decode('utf-8')
            ndata += 1
            if ndata % DINKUMCHUNKSIZE == 0:
                # need to allocate more data!
                data = np.concatenate(
                    (data, np.NaN + np.zeros((DINKUMCHUNKSIZE, nsensors))),
                    axis=0)
        elif d == 'X':
            # End of file cycle tag. We made it through.
            # throw out pre-allocated data we didn't use...
            data = data[:ndata]
            break
        else:
            raise ValueError(('Parsing failed at {}. ',
                'Got {} expected d or X').format(binaryData.bytepos, d))

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

    # Calculate flight times for science data
    uniqueSciTimes, uniqueSciTimeIndices = np.unique(
            np.array(fdata['sci_m_present_time']), return_index=True)
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
            sdata['m_present_time_sci'] = np.NaN * sdata['sci_m_present_time']

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

def _make_dinkumcache(filelist, cachedir):
    """
    Helper function to setup the cache of sensor names based on the crc
    number in the header of the first file in filelist
    """
    for filen in filelist:
        try:
            # keep trying files until we find one that makes the cache...
            meta, bindatafilepos = dbd_get_meta(filen, cachedir)
            return True
        except:
            pass
    return False


def datameta_to_nc(data, meta, outdir=None, name=None, check_exists=False,
                deployment_ind=0):
    """
    Convert a raw dinkum data and meta dict to a netcdf fileself.

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

    ds['time'] = (('_ind'), time, {'units': 'seconds since 1970-01-01T00:00:00Z'})
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
    _log.debug(f'{ds}, {ds.time}, {len(ind)}')
    ds = ds.isel(_ind=ind)
    ds['_ind'] = np.arange(len(ds.time)) + deployment_ind
    if len(ds['_ind'].values) > 1:
        deployment_ind = ds['_ind'].values[-1]

    _log.info(f'Writing! {deployment_ind}')
    ds.to_netcdf(outname, format='NETCDF4')
    _log.info(f'Wrote:, {outname}')
    return ds, deployment_ind

def merge_rawnc(indir, outdir, deploymentyaml, incremental=False,
                scisuffix='EBD', glidersuffix='DBD'):

    scisuffix = scisuffix.lower()
    glidersuffix = glidersuffix.lower()
    with open(deploymentyaml) as fin:
        deployment = yaml.safe_load(fin)
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

    for num in range(1, 500):
        d = indir + f'/{num:04d}*.' + glidersuffix + '.nc'
        _log.info('Check %s', d)
        fin = glob.glob(d)
        if fin:
            d = indir + f'/{num:04d}*.' + glidersuffix + '.nc'
            if glob.glob(d):
                outnebd = outdir + '/' + id + f'-{num:04d}-rawdbd.nc'
                _log.info('Opening *.dbd.nc multi-file dataset')
                _log.info(outnebd)

                with xr.open_mfdataset(d, decode_times=False, lock=False) as ds:
                    ds = ds.sortby('time')
                    ds['_ind'] = np.arange(len(ds.time))
                    ds.to_netcdf(outnebd, 'w')

        d = indir + f'/{num:04d}*.' + scisuffix + '.nc'
        fin = glob.glob(d)
        if fin:
            d = indir + f'/{num:04d}*.' + scisuffix + '.nc'
            if glob.glob(d):
                outndbd = outdir + '/' + id + f'-{num:04d}-rawebd.nc'
                _log.info('Opening *.ebd.nc multi-file dataset')
                with xr.open_mfdataset(d, decode_times=False,  lock=False) as ds:
                    ds = ds.sortby('time')
                    ds['_ind'] = np.arange(len(ds.time))
                    ds.to_netcdf(outndbd, 'w')


def merge_rawncBrutal(indir, outdir, deploymentyaml, incremental=False,
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

    incremental : bool
        Only add new files....
    """

    with open(deploymentyaml) as fin:
        deployment = yaml.safe_load(fin)
    metadata = deployment['metadata']
    id = metadata['glider_name'] + metadata['glider_serial']
    outnebd = outdir + '/' + id + '-rawebd.nc'
    outndbd = outdir + '/' + id + '-rawdbd.nc'

    _log.info('Opening *.ebd.nc multi-file dataset')
    scifiles = sorted(glob.glob(indir + '/*.' + scisuffix + '.nc'))
    for fn in scifiles:
        _log.debug(fn)
    glifiles = sorted(glob.glob(indir + '/*.' + glidersuffix + '.nc'))

    if len(scifiles) > 1:
        # this is complicated because merging is very slow compared to
        # mfd dataset.  So merge 10 at a time, save the intermeddairies,
        # and then re-open.  We have to do this file by file because some
        # of the files are bad...
        ds = None
        num = 0
        for nn, name in enumerate(scifiles):
            _log.debug('Hello')
            _log.info(f'merging {scifiles[nn]}')
            with xr.open_dataset(name, decode_times=True) as ds2:
                if 1:
                    if ds is not None:
                        ds = ds.merge(ds2)
                    else:
                        ds = ds2.copy()
                    num = num + 1
                else:
                    _log.info(f'Failed to merge {name}')
            if num % 10 == 0 or nn == len(scifiles) - 1:
                ds.to_netcdf(indir+f'TEMP{num:04d}.nc', 'w')
                ds = None
        with xr.open_mfdataset(indir+f'TEMP*.nc', decode_times=False, lock=False) as ds:
            dsnew = ds.sortby(ds.time)
            dsnew.to_netcdf(outnebd, 'w')
            _log.info('Wrote ' + outnebd)

    if len(glifiles) > 1:
        ds = None
        num = 0
        for nn, name in enumerate(glifiles):
            _log.info(f'merging {glifiles[nn]}')
            with xr.open_dataset(name, decode_times=False) as ds2:
                try:
                    if ds is not None:
                        ds = ds.merge(ds2)
                    else:
                        ds = ds2.copy()
                    num = num + 1
                except:
                    _log.info(f'Failed to merge {name}')
            if num % 10 == 0 or nn == len(scifiles) - 1:
                ds.to_netcdf(indir+f'TEMPG{num:04d}.nc', 'w')
                ds = None
        with xr.open_mfdataset(indir+f'TEMPG*.nc', decode_times=False, lock=False) as ds:
            dsnew = ds.sortby(ds.time)
            dsnew.to_netcdf(outndbd, 'w')
            _log.info('Wrote ' + outndbd)

    return


#only the merge jpnote

## jpnote: global variables ##
# - preprocess overwrites variables -
# - need global list to append multiple netcdf file attributes ## 

globlist = []
latmin = []
latmax = []
lonmin = []
lonmax = []

def preprocess(ds):

    global globlist
    global latmin
    global latmax
    global lonmin
    global lonmax

  #  dsnew = ds.copy()
# get time coverage start and time coverage end from each necdf file 
    startmin = np.datetime64(ds.attrs['time_coverage_start']) 
    endmax = np.datetime64(ds.attrs['time_coverage_end'])
# append each to global list for later min max extraction 
    globlist.append(startmin)
    globlist.append(endmax)
# get geospatial min and max from each necdf file 
    latminv = ds.attrs['geospatial_lat_min'] 
    latmaxv = ds.attrs['geospatial_lat_max']
    lonminv = ds.attrs['geospatial_lon_min']
    lonmaxv = ds.attrs['geospatial_lon_max']
# append each to global list for later min max extraction 
    latmin.append(latminv)
    latmax.append(latmaxv)
    lonmin.append(lonminv)
    lonmax.append(lonmaxv)

    return ds

##jpnote undo this comment after testing 

def raw_to_L0timeseries(indir, outdir, deploymentyaml, *,
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

    #jpnote: for max min extraction for attributes of merged file 
    global globlist
    global latmax
    global latmin 
    global lonmin
    global lonmax 

    with open(deploymentyaml) as fin:
        deployment = yaml.safe_load(fin)
    metadata = deployment['metadata']
    ncvar = deployment['netcdf_variables']
    device_data = deployment['glider_devices']
    thenames = list(ncvar.keys())
    thenames.remove('time')

    id = metadata['glider_name'] + metadata['glider_serial']

    id0 = None
    prev_profile = 0
    for mnum in range(0, 500):
        if 1:
            ebdn = indir + '/' + id + f'-{mnum:04d}-rawebd.nc'
            dbdn = indir + '/' + id + f'-{mnum:04d}-rawdbd.nc'
            _log.debug(f'{ebdn}, {dbdn}')
            if os.path.exists(ebdn) and os.path.exists(dbdn):
                _log.info('Opening:', ebdn, dbdn)
                ebd = xr.open_dataset(ebdn, decode_times=False)
                dbd = xr.open_dataset(dbdn, decode_times=False)
                _log.debug(f'DBD, {dbd}, {dbd.m_depth}')
                if len(ebd.time) > 2:
                    # build a new data set based on info in `deployment.`
                    # We will use ebd.m_present_time as the interpolant if the
                    # variabel is in dbd.

                    ds = xr.Dataset()
                    attr = {}
                    name = 'time'
                    for atts in ncvar[name].keys():
                        if atts != 'coordinates':
                            attr[atts] = ncvar[name][atts]
                    ds[name] = (('time'), ebd[name].values, attr)

                    for name in thenames:
                        _log.info('working on %s', name)
                        if not('method' in ncvar[name].keys()):
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
                        #        val[val==0] = np.NaN
                                val = convert(val)
                            else:
                                _log.debug('DBD sensorname %s', sensorname)
                                val = convert(dbd[sensorname])
                                val = _dbd2ebd(dbd, ds, val)
                                ncvar['method'] = 'linear fill'
                            # make the attributes:
                            ncvar[name].pop('coordinates', None)
                            attrs = ncvar[name]
                            attrs = utils.fill_required_attrs(attrs)
                            ds[name] = (('time'), val.data, attrs)

                    _log.debug(f'HERE, {ds}')
                    _log.debug(f'HERE, {ds.pressure[0:100]}')
                    # some derived variables:
                    # trim bad times...
                    #ds = ds.sel(time=slice(1e8, None))

                    ds = utils.get_glider_depth(ds)

                    ds = utils.get_distance_over_ground(ds)

                    # ds = utils.get_profiles(ds)
                    # ds['profile_index'] = ds.profile_index + prev_profile

                    #ind = np.where(np.isfinite(ds.profile_index))[0]
                    #prev_profile = ds.profile_index.values[ind][-1]
                    ds = utils.get_derived_eos_raw(ds)
                    ds = ds.assign_coords(longitude=ds.longitude)
                    ds = ds.assign_coords(latitude=ds.latitude)
                    ds = ds.assign_coords(depth=ds.depth)

                    #ds = ds._get_distance_over_ground(ds)
                    ds['time'] = ds.time.values.astype('timedelta64[s]') + np.datetime64('1970-01-01T00:00:00')
                    ds = utils.fill_metadata(ds, deployment['metadata'], device_data)
                    try:
                        os.mkdir(outdir)
                    except:
                        pass
                    outname = (outdir + '/' + ds.attrs['deployment_name'] +
                                f'-M{mnum:04d}_L0.nc')
                    _log.info('writing %s', outname)
                    ds.to_netcdf(outname, 'w', encoding={'time': {'units': 'seconds since 1970-01-01T00:00:00Z'}})
                    if id0 is None:
                        id0 = ds.attrs['deployment_name']

# now merge:
    #jpnote: added preprocess function to xr.open_mfdataset
    with xr.open_mfdataset(outdir + '/' + id + '*-M*_L0.nc', lock=False,preprocess=preprocess) as ds:
        _log.debug(ds.attrs)

        # put the real start and end times:
        start = ds['time'].values[0]
        end = ds['time'].values[-1]

        ds.attrs['deployment_start'] = str(start)
        ds.attrs['deployment_end'] = str(end)

## jpnote: beginning of jp addition ## 
        time_coverage_start_jp = str(min(globlist))
        time_coverage_end_jp = str(max(globlist))
# add calcs to the dataset 
        ds.attrs['time_coverage_start'] = time_coverage_start_jp 
        ds.attrs['time_coverage_end'] = time_coverage_end_jp
# remove 0.0 from list, so that proper min,max can be extracted
        latmin = [i for i in latmin if i != 0]
        latmax = [i for i in latmax if i != 0]
        lonmin = [i for i in lonmin if i != 0]
        lonmax = [i for i in lonmax if i != 0]
# add calcs to dataset 
        ds.attrs['geospatial_lat_min'] = min(latmin)
        ds.attrs['geospatial_lat_max'] = max(latmax)
        ds.attrs['geospatial_lon_min'] = min(lonmin)
        ds.attrs['geospatial_lon_max'] = max(lonmax)

##end of jp addition 
        _log.debug(ds.depth.values[:100])
        _log.debug(ds.depth.values[2000:2100])
        ds = utils.get_profiles_new(ds,
                filt_time=profile_filt_time, profile_min_time=profile_min_time)
        _log.debug(ds.depth.values[:100])
        _log.debug(ds.depth.values[2000:2100])

        outname = outdir + '/' + id0 + '.nc'
        _log.info(outname)
        ds.to_netcdf(outname)
    return outname

# alias:
raw_to_L1timeseries = raw_to_L0timeseries

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
    with  xr.open_dataset(inname) as ds:
        ds = utils.get_profiles_new(ds,
                filt_time=profile_filt_time, profile_min_time=profile_min_time)
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
        vout[goodt] = np.interp(ds.time[goodt].values,
            dbd.m_present_time.values[good], val[good].values)
    return vout
