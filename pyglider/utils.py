import collections
import seawater
import xarray as xr
import numpy as np
from scipy.signal import argrelextrema
# import webcolors
import logging

_log = logging.getLogger(__name__)


def get_distance_over_ground(ds):
    good = ~np.isnan(ds.latitude + ds.longitude)
    dist, cog = seawater.dist(ds.latitude[good], ds.longitude[good])
    dist = np.roll(np.append(dist, 0), 1)
    dist = np.cumsum(dist)
    attr = {'long_name': 'distance over ground flown since mission start',
     'method': 'get_distance_over_ground',
     'units': 'km',
     'sources': 'latitude longitude'}
    ds['distance_over_ground'] = (('time'), dist, attr)
    return ds

def get_glider_depth(ds):

    good = np.where(~np.isnan(ds.pressure))[0]
    ds['depth'] = ds.pressure * 0.
    ds['depth'].values = seawater.eos80.dpth(ds.pressure.values,
            ds.latitude.mean().values)
    # now we really want to know where it is, so interpolate:
    if len(good) > 0:
        ds['depth'].values = np.interp(np.arange(len(ds.depth)),
                                good, ds['depth'].values[good])

    attr = {'source': 'pressure', 'long_name': 'glider depth',
            'standard_name': 'depth', 'units': 'm',
            'comment': 'from science pressure and interpolated',
            'observation_type': 'calulated',
            'accuracy': '1', 'precision': '2', 'resolution': '0.02',
            'valid_min': '0', 'valid_max': '2000',
            'reference_datum': 'surface', 'positive': 'down'}
    ds['depth'].attrs = attr
    return ds


def get_profiles(ds, min_dp = 10.0, inversion=3., filt_length=7,
                 min_nsamples=14):
    """
    make two variables: profile_direction and profile_index; this veersion
    is good for lots of data.  Less good for sparse data
    """
    profile = ds.pressure.values * 0
    direction = ds.pressure.values * 0
    pronum = 1
    lastpronum = 0

    good = np.where(~np.isnan(ds.pressure))[0]
    p = np.convolve(ds.pressure.values[good],
                    np.ones(filt_length) / filt_length, 'same')
    dpall = np.diff(p)
    inflect = np.where(dpall[:-1] * dpall[1:] < 0)[0]
    for n, i in enumerate(inflect[:-1]):
        nprofile = inflect[n+1] - inflect[n]
        inds = np.arange(good[inflect[n]], good[inflect[n+1]]+1) + 1
        dp = np.diff(ds.pressure[inds[[-1, 0]]])
        if ((nprofile >= min_nsamples) and (np.abs(dp) > 10)):
            print('Good')
            direction[inds] = np.sign(dp)
            profile[inds] = pronum
            lastpronum = pronum
            pronum += 1
        else:
            profile[good[inflect[n]]:good[inflect[n+1]]] = lastpronum + 0.5

    attrs = collections.OrderedDict([('long_name', 'profile index'),
             ('units', '1'),
             ('comment',
              'N = inside profile N, N + 0.5 = between profiles N and N + 1'),
             ('sources', 'time pressure'),
             ('method', 'get_profiles'),
             ('min_dp', min_dp),
             ('inversion', inversion),
             ('filt_length', filt_length),
             ('min_nsamples', min_nsamples)])
    ds['profile_index'] = (('time'), profile, attrs)


    attrs = collections.OrderedDict([('long_name', 'glider vertical speed direction'),
             ('units', '1'),
             ('comment',
              '-1 = ascending, 0 = inflecting or stalled, 1 = descending'),
             ('sources', 'time pressure'),
             ('method', 'get_profiles')])
    ds['profile_direction'] = (('time'), direction, attrs)
    return ds


def get_profiles_new(ds, min_dp = 10.0, inversion=3., filt_time=100,
                 profile_min_time=300):
    """
    make two variables: profile_direction and profile_index; this version
    is good for lots of data.  Less good for sparse data

    filt_time is in seconds
    min_nsamples is in seconds
    """
    profile = ds.pressure.values * 0
    direction = ds.pressure.values * 0
    pronum = 1
    lastpronum = 0

    good = np.where(np.isfinite(ds.pressure))[0]
    dt = float(np.median(np.diff(ds.time.values[good[:200000]])))
    print('dt', dt)
    filt_length = int(filt_time /  dt)

    min_nsamples = int(profile_min_time / dt)
    _log.info('Filt Len  %d, dt %f, min_n %d', filt_length, dt, min_nsamples)

    p = np.convolve(ds.pressure.values[good],
                    np.ones(filt_length) / filt_length, 'same')
    print('filt', filt_length)
    decim = int(filt_length / 3)
    if decim < 2:
        decim = 2
    # why?  because argrelextrema doesn't like repeated values, so smooth
    # then decimate to get fewer values:
    pp = p[::decim]
    maxs = argrelextrema(pp, np.greater)[0]
    mins = argrelextrema(pp, np.less)[0]
    mins = good[mins * decim]
    maxs = good[maxs * decim]
    if mins[0] > maxs[0]:
        mins = np.concatenate(([0], mins))
    if mins[-1] < maxs[-1]:
        mins = np.concatenate((mins, good[[-1]]))

    _log.info(f'mins: {len(mins)} {mins} , maxs: {len(maxs)} {maxs}')

    pronum = 0
    p = ds.pressure
    nmin = 0
    nmax=0
    while (nmin < len(mins)) and (nmax < len(maxs)):
        if 1:
            nmax = np.where(maxs>mins[nmin])[0]
            if len(nmax) >= 1:
                nmax = nmax[0]
            else:
                break
            print(nmax)
            ins = range(int(mins[nmin]), int(maxs[nmax]+1))
            print(pronum, ins, len(p), mins[nmin], maxs[nmax])
            print('Down', ins, p[ins[0]].values,p[ins[-1]].values)
            if (len(ins) > min_nsamples and np.nanmax(p[ins]) - np.nanmin(p[ins]) > min_dp):
                profile[ins] = pronum
                direction[ins] = +1
                pronum += 1
            nmin = np.where(mins>maxs[nmax])[0]
            if len(nmin) >= 1:
                nmin = nmin[0]
            else:
                break
            ins = range(maxs[nmax], mins[nmin])
            print(pronum, ins, len(p), mins[nmin], maxs[nmax])
            print('Up', ins, p[ins[0]].values, p[ins[-1]].values)
            if (len(ins) > min_nsamples and np.nanmax(p[ins]) - np.nanmin(p[ins]) > min_dp):
                # up
                profile[ins] = pronum
                direction[ins] = -1
                pronum += 1
        else:
            print('Failed?')

    print('Doing this...')
    attrs = collections.OrderedDict([('long_name', 'profile index'),
             ('units', '1'),
             ('comment',
              'N = inside profile N, N + 0.5 = between profiles N and N + 1'),
             ('sources', 'time pressure'),
             ('method', 'get_profiles'),
             ('min_dp', min_dp),
             ('inversion', inversion),
             ('filt_length', filt_length),
             ('min_nsamples', min_nsamples)])
    ds['profile_index'] = (('time'), profile, attrs)


    attrs = collections.OrderedDict([('long_name', 'glider vertical speed direction'),
             ('units', '1'),
             ('comment',
              '-1 = ascending, 0 = inflecting or stalled, 1 = descending'),
             ('sources', 'time pressure'),
             ('method', 'get_profiles')])
    ds['profile_direction'] = (('time'), direction, attrs)
    return ds


def get_derived_eos_raw(ds):
    # GPCTD and slocum ctd require a scale factor of 10 for conductivity. Legato does not
    if 'GPCTD' in ds.conductivity.source or 'sci_water_cond' in ds.conductivity.source:
        ds['conductivity'] = ds['conductivity'] * 10
    r = ds.conductivity / seawater.constants.c3515
    ds['salinity'] = (('time'),
                      seawater.eos80.salt(r, ds.temperature, ds.pressure))
    attrs = collections.OrderedDict([('long_name', 'water salinity'),
             ('standard_name', 'sea_water_practical_salinity'),
             ('units', '1e-3'),
             ('comment', 'raw, uncorrected salinity'),
             ('sources', 'conductivity temperature pressure'),
             ('method', 'get_derived_eos_raw'),
             ('observation_type', 'calulated'),
             ('instrument', 'instrument_ctd'),
             ('valid_max', '40.0'),
             ('valid_min', '0.0'),
             ('accuracy', '0.01'),
             ('precision', '0.01'),
             ('resolution', '0.001')
             ])
    attrs = fill_required_attrs(attrs)
    ds['salinity'].attrs = attrs

    ds['potential_density'] = (('time'), seawater.eos80.pden(
            ds.salinity, ds.temperature, ds.pressure, pr=0))
    attrs = collections.OrderedDict([('long_name', 'water potential density'),
             ('standard_name', 'sea_water_potential_density'),
             ('units', 'kg m-3'),
             ('comment', 'raw, uncorrected salinity'),
             ('sources', 'salinity temperature pressure'),
             ('method', 'get_derived_eos_raw'),
             ('observation_type', 'calulated'),
             ('instrument', 'instrument_ctd'),
              ('accuracy', '0.01'),
              ('precision', '0.01'),
              ('resolution', '0.001')
             ])
    attrs = fill_required_attrs(attrs)
    ds['potential_density'].attrs = attrs

    ds['density'] = (('time'), seawater.eos80.dens(
            ds.salinity, ds.temperature, ds.pressure))
    attrs = collections.OrderedDict([('long_name', 'Density'),
             ('standard_name', 'sea_water_density'),
             ('units', 'kg m-3'),
             ('comment', 'raw, uncorrected salinity'),
             ('observation_type', 'calulated'),
             ('sources', 'salinity temperature pressure'),
             ('instrument', 'instrument_ctd'),
             ('method', 'get_derived_eos_raw'),
              ('valid_min', '1000.0'),
              ('valid_max', '1040.0'),
               ('accuracy', '0.01'),
               ('precision', '0.01'),
               ('resolution', '0.001')
             ])
    attrs = fill_required_attrs(attrs)
    ds['density'].attrs = attrs

    ds['potential_temperature'] = (('time'), seawater.eos80.ptmp(
            ds.salinity, ds.temperature, ds.pressure, pr=0))
    attrs = collections.OrderedDict([('long_name', 'water potential temperature'),
             ('standard_name', 'sea_water_potential_temperature'),
             ('units', 'Celsius'),
             ('comment', 'raw, uncorrected salinity'),
             ('sources', 'salinity temperature pressure'),
             ('observation_type', 'calulated'),
             ('method', 'get_derived_eos_raw'),
              ('instrument', 'instrument_ctd'),
           ('accuracy', '0.002'),
           ('precision', '0.001'),
           ('resolution', '0.0001')
             ])
    attrs = fill_required_attrs(attrs)
    ds['potential_temperature'].attrs = attrs

    return ds


def time_to_datetime64(time):
    """
    Pass in a glider undecoded time (seconds since 1970-01-01), and
    get a np.datetime64[s] object back.
    """
    return (time.astype('timedelta64[s]') +
            np.datetime64('1970-01-01'))


def fill_required_attrs(attrs):
    required = {
        'comment': " ",
        'accuracy': " ",
        'precision': " ",
        'platform':  "platform",
        'resolution': " ",
        'ancillary_variables': " "}
    for k in required.keys():
        if not k in attrs.keys():
            attrs[k] = required[k]
    return attrs


def get_file_id(ds):

    print(ds.time)
    if not ds.time.dtype=='datetime64[ns]':
        dt = time_to_datetime64(ds.time.values[0])
    else:
        dt = ds.time.values[0].astype('datetime64[s]')
    print('dt', dt)
    id = (ds.attrs['glider_name'] + ds.attrs['glider_serial'] + '-' +
                      dt.item().strftime('%Y%m%dT%H%M'))
    return id


def fill_metadata(ds, metadata):

    good = ~np.isnan(ds.latitude.values + ds.longitude.values)
    ds.attrs['geospatial_lat_max'] = np.max(ds.latitude.values[good])
    ds.attrs['geospatial_lat_min'] = np.min(ds.latitude.values[good])
    ds.attrs['geospatial_lon_max'] = np.max(ds.longitude.values[good])
    ds.attrs['geospatial_lon_min'] = np.min(ds.longitude.values[good])
    ds.attrs['geospatial_lat_units'] = 'degrees_north'
    ds.attrs['geospatial_lon_units'] = 'degrees_east'
    ds.attrs['netcdf_version'] = '4.0'  # TODO get this somehow...
    ds.attrs['history'] = 'CPROOF glider toolbox version: pre-tag'
    for k, v in metadata.items():
        ds.attrs[k] = v
    ds.attrs['featureType'] = 'timeseries'
    ds.attrs['cdm_data_type'] = 'Trajectory'
    ds.attrs['Conventions'] = 'CF-1.6'
    ds.attrs['date_created'] = str(np.datetime64('now')) + 'Z'
    ds.attrs['date_issued'] = str(np.datetime64('now')) + 'Z'
    ds.attrs['date_modified'] = " "
    ds.attrs['id'] = get_file_id(ds)
    ds.attrs['deployment_name'] = metadata['deployment_name']
    ds.attrs['title'] = ds.attrs['id']

    dt = time_to_datetime64(ds.time.values)
    ds.attrs['time_coverage_start'] = '%s' % dt[0]
    ds.attrs['time_coverage_end'] = '%s' % dt[-1]

    ds.attrs['processing_level'] = 'Level 0 (L0) processed data timeseries; no corrections or data screening'

    ds.attrs = collections.OrderedDict(sorted(ds.attrs.items()))

    return ds

def _zero_screen(val):
    val[val==0] = np.NaN
    return val


def nmea2deg(nmea):
    deg = (np.fix(nmea / 100) +
           np.sign(nmea) * np.remainder(np.abs(nmea), 100) / 60)
    return deg

def bar2dbar(val):
    return val * 10.0

def cbar2dbar(val):
    return val / 10.0


def _passthrough(val):
    return val


def gappy_fill_vertical(data):
    """
    Fill vertical gaps from the first to last bin with data in them.
    Applied column-wise.

    data = gappy_fill_vertical(data)
    """
    m, n = np.shape(data)
    for j in range(n):
        ind = np.where(~np.isnan(data[:, j]))[0]
        if (len(ind) > 0 and len(ind) < (ind[-1] - ind[0])
                and len(ind) > (ind[-1] - ind[0]) * 0.05):
            int = np.arange(ind[0], ind[-1])
            data[:, j][ind[0]:ind[-1]] = np.interp(int, ind, data[ind, j])
    return data

def parse_gliderxml_pos(fname):
    """
    returns lon, lat, timestr
    """

    xmln = fname
    from bs4 import BeautifulSoup
    with open(xmln, 'r') as fin:
        y=BeautifulSoup(fin, features="xml")
        time = None
        lon = None
        lat = None
        for a in y.find_all('valid_location'):
            try:
                dtime = np.datetime64(a.time.text)
                if dtime > np.datetime64('2019-07-19'):
                    #print(dtime, float(a.lat.text), float(a.lon.text))
                    #print(np.array(dtime))
                    if time is not None:
                        time = np.append(time, dtime)
                        lat = np.append(lat, float(a.lat.text))
                        lon = np.append(lon, float(a.lon.text))
                    else:
                        time = np.array(dtime)
                        lat = np.array(float(a.lat.text))
                        lon = np.array(float(a.lon.text))
            except:
                pass
        lon = nmea2deg(lon)
        lat = nmea2deg(lat)
    return lon, lat, time

def parse_gliderxml_surfacedatetime(fname):
    """
    returns lon, lat, timestr
    """

    xmln = fname
    from bs4 import BeautifulSoup
    with open(xmln, 'r') as fin:
        y=BeautifulSoup(fin, features="xml")
        time = None
        for a in y.find_all('report'):
            print(a)
            if a.text is not None:
                try:
                    time = np.append(time, np.datetime64(a.text))
                except:
                    pass
    return time

def get_html_non_blue(num=None):
    """
    """

    colors = list(webcolors.css3_names_to_hex.keys())
    colors0 = colors.copy()
    for c in colors0:
        if (c.lower().find('blue') >= 0 or
                c.lower().find('cyan') >= 0 or
                c.lower().find('turquoise') >= 0 or
                c.lower().find('aqua')>= 0 or
                c.lower().find('navy') >= 0 or
                c.lower().find('teal')>= 0 or
                c.lower().find('teal')>= 0 or
                c.lower().find('gray')>= 0 or
                c.lower().find('grey')>= 0  or
                c.lower().find('white')>= 0  or
                c.lower().find('sea')>= 0):
            colors.remove(c)
    if num is None:
        num = np.random.randint(0, len(colors))
    cname = colors[num]
    ints = webcolors.name_to_rgb(cname)
    return ints
