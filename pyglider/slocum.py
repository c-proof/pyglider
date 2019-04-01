# -*- coding: utf-8 -*-
from pyglider import bitstring
import datetime
import itertools
import logging
from math import floor, fmod
import numpy as np
import os
import re
import time
import xarray as xr


_log = logging.getLogger(__name__)


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
    fname = cachedir + '/' + meta['sensor_list_crc'].upper() + '.CAC'
    with open(fname, 'rb') as dfh:
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


def dbd_get_meta(filename, cachedir=None):
    """
    Get metadata from a dinkum binary file.

    Parameters
    ----------

    filename : str
        filename of the dinkum binary file (i.e. *.dbd, *.ebd)

    cachedir : str or None
        Directory where the cached CAC sensor lists are kept.  These
        lists are often in directories like ``../Main_board/STATE/CACHE/``.
        These should be copied somewhere locally.  If None, defaults to and
        creates ``./.dinkumcache`` in the local directory.

    Returns
    -------
    meta : dict
        Dictionary of the meta data for this dinkum binary file.

    """
    if cachedir is None:
        cachedir = './.dinkumcache/'

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


def dbd_to_dict(dinkum_file, cachedir=None, keys=None):
    """
    Translate a dinkum binary file to a dictionary of data and meta values.

    Parameters
    ----------
    dinkum_file : dbd file name (full path)
        These are the raw data from the glider, either offloaded from a card
        or from the dockserver.

    cachedir : str or None
        Directory where the cached CAC sensor lists are kept.  These
        lists are often in directories like ``../Main_board/STATE/CACHE/``.
        These should be copied somewhere locally.  If None, defaults to and
        creates ``./.dinkumcache`` in the local directory.

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
    if cachedir is None:
        cachedir = './.dinkumcache/'

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
            if code == '01':  # Same value as before.
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

def _make_dinkumcache(filelist, cachedir=None):
    """
    Helper function to setup the cache of sensor names based on the crc
    number in the header of the first file in filelist
    """
    if cachedir is None:
        cachedir = './.dinkumcache/'
    for filen in filelist:
        try:
            # keep trying files until we find one that makes the cache...
            meta, bindatafilepos = dbd_get_meta(filen, cachedir)
            return True
        except:
            pass
    return False


def datameta_to_nc(data, meta, outdir=None, name=None, check_exists=False):
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

    ds['time'] = (('time'), time, {'units': 'seconds since 1970-01-01T00:00:00Z'})
    for key in data.keys():
        ds[key] = (('time'), data[key])
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

    ds.to_netcdf(outname, format='NETCDF4')
    return ds


def _mergeMultiple(filelist, cachedir=None, keys=None):
    """
    Not currently used....

    Merge the dinkum ebd/dbd files in filelist into a dictionary with each
    time series as a key name.

    These should all be the same type of file, ie Science files (*.ebd) or
    Main_board files (*.dbd).  Mixing them won't work.

    Parameters:
    -----------
    filelist : list of files to merge
        i.e. ``filelist = glob.glob('datadirectory/Science/*.ebd')``

    cachedir : directory to read/store sensor list caches
        If not specified, the cache is assumed to be kept locally
        in ``.dinkumcache/``.

    keys : list of strings or a str
        If a list of strings just return dictionary with those strings.
        If a single string use `dinkum.parse_filter_file` to get the
        keys

    Returns:
    --------
    outdata : dictionary of data

    outmeta : list of meta data from each file merged.

    Notes:
    ------

    Most files will not have a list of sensors;  These can be read from
    a cached file containing the list, or the files in *filelist* are
    searched until a list of sensors is found.

    outdata = dict()
    outmeta = []

    if isinstance(keys, str):
        keys = parse_filter_file(keys)

    for n in range(len(filelist)):
        f = filelist[n]
        fdata, fmeta = dbd_to_dict(f, keys=keys)
        for key, dd in fdata.items():
            if key in outdata:
                outdata[key] = np.concatenate((outdata[key], dd))
            else:
                outdata[key] = dd
        outmeta.append(fmeta)

    return outdata, outmeta
    """
    pass


def _webb_to_decdeg(webb_latlon):
    """
    Not currently used.  Maybe useful later...
    if webb_latlon is None or webb_latlon == 'NaN' or webb_latlon == '69696969.0':
        return None
    deg = floor(abs(webb_latlon) / 100)
    minute = fmod(abs(webb_latlon), 100) / 60
    if webb_latlon > 0:
        return deg + minute
    else:
        return -1 * (abs(deg) + minute)
    """
    pass
