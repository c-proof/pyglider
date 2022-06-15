import collections
import xarray as xr
import numpy as np
from scipy.signal import argrelextrema
import gsw
# import webcolors
import logging

_log = logging.getLogger(__name__)


def get_distance_over_ground(ds):
    good = ~np.isnan(ds.latitude + ds.longitude)
    dist = gsw.distance(ds.longitude[good].values, ds.latitude[good].values)/1000
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
    ds['depth'].values = -gsw.z_from_p(ds.pressure.values,
            ds.latitude.values)
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
            _log.debug('Good')
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
    #dt = float(np.median(np.diff(ds.time.values[good[:200000]])))
    dt = float(np.median(np.diff(ds.time.values[good[:200000]]).astype(np.float64)) * 1e-9)
    _log.info(f'dt, {dt}')
    filt_length = int(filt_time /  dt)

    min_nsamples = int(profile_min_time / dt)
    _log.info('Filt Len  %d, dt %f, min_n %d', filt_length, dt, min_nsamples)

    p = np.convolve(ds.pressure.values[good],
                    np.ones(filt_length) / filt_length, 'same')
    _log.info('ffilt, {filt_length}')
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

    _log.debug(f'mins: {len(mins)} {mins} , maxs: {len(maxs)} {maxs}')

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
            _log.debug(nmax)
            ins = range(int(mins[nmin]), int(maxs[nmax]+1))
            _log.debug(f'{pronum}, {ins}, {len(p)}, {mins[nmin]}, {maxs[nmax]}')
            _log.debug(f'Down, {ins}, {p[ins[0]].values},{p[ins[-1]].values}')
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
            _log.debug(f'{pronum}, {ins}, {len(p)}, {mins[nmin]}, {maxs[nmax]}')
            _log.debug(f'Up, {ins}, {p[ins[0]].values}, {p[ins[-1]].values}')
            if (len(ins) > min_nsamples and np.nanmax(p[ins]) - np.nanmin(p[ins]) > min_dp):
                # up
                profile[ins] = pronum
                direction[ins] = -1
                pronum += 1
        else:
            _log.debug('Failed?')

    _log.debug('Doing this...')
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
    if 'S m' in ds.conductivity.units:
        r = 10 * ds.conductivity
    elif 'mS cm' in ds.conductivity.units:
        r = ds.conductivity
    else:
        raise ValueError("Could not parse conductivity units in yaml. Expected 'S m-1' or 'mS cm-1'. "
                         "Check yaml entry netcdf_variables: conductivity: units")
    ds['salinity'] = (('time'),
                      gsw.conversions.SP_from_C(r, ds.temperature, ds.pressure).values)
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
    sa = gsw.SA_from_SP(ds['salinity'], ds['pressure'], ds['longitude'], ds['latitude'])
    ct = gsw.CT_from_t(sa, ds['temperature'], ds['pressure'])
    ds['potential_density'] = (('time'), 1000 + gsw.density.sigma0(sa, ct).values)
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

    ds['density'] = (('time'), gsw.density.rho(
            ds.salinity, ds.temperature, ds.pressure).values)
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
    ds['potential_temperature'] = (('time'), gsw.conversions.pt0_from_t(
            ds.salinity, ds.temperature, ds.pressure).values)
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

    _log.debug(ds.time)
    if not ds.time.dtype=='datetime64[ns]':
        dt = time_to_datetime64(ds.time.values[0])
    else:
        dt = ds.time.values[0].astype('datetime64[s]')
    _log.debug(f'dt, {dt}')
    id = (ds.attrs['glider_name'] + ds.attrs['glider_serial'] + '-' +
                      dt.item().strftime('%Y%m%dT%H%M'))
    return id


def fill_metadata(ds, metadata, sensor_data):

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

    dt = ds.time.values
    ds.attrs['time_coverage_start'] = '%s' % dt[0]
    ds.attrs['time_coverage_end'] = '%s' % dt[-1]

    ds.attrs['processing_level'] = 'Level 0 (L0) processed data timeseries; no corrections or data screening'

    for k, v in sensor_data.items():
        ds.attrs[k] = str(v)
    ds.attrs = collections.OrderedDict(sorted(ds.attrs.items()))

    return ds

def _zero_screen(val):
    val[val==0] = np.NaN
    return val


def nmea2deg(nmea):
    deg = (np.fix(nmea / 100) +
           np.sign(nmea) * np.remainder(np.abs(nmea), 100) / 60)
    return deg


def oxygen_concentration_correction(data, ncvar):
    oxy_yaml = ncvar['oxygen_concentration']
    if 'reference_salinity' not in oxy_yaml.keys():
        _log.warning('No reference_salinity found in oxygen deployment yaml. Assuming reference salinity of 0 psu')
        ref_sal = 0
    else:
        ref_sal = float(oxy_yaml['reference_salinity'])
    _log.info(f'Correcting oxygen using reference salinity {ref_sal} PSU')
    ds_oxy = data.oxygen_concentration[~np.isnan(data.oxygen_concentration)]
    # Match the nearest temperature and salinity values from their timestamps
    ds_temp = data.potential_temperature[~np.isnan(data.potential_temperature)].reindex(time=ds_oxy.time, method="nearest")
    ds_sal = data.salinity[~np.isnan(data.salinity)].reindex(time=ds_oxy.time, method="nearest")
    o2_sol = gsw.O2sol_SP_pt(ds_sal, ds_temp)
    o2_sat = ds_oxy / gsw.O2sol_SP_pt(ds_sal*0 + ref_sal, ds_temp)
    data['oxygen_concentration'].values[~np.isnan(data.oxygen_concentration)] = o2_sat * o2_sol
    data['oxygen_concentration'].attrs['oxygen_concentration_QC:RTQC_methodology'] =\
        f'oxygen concentration corrected for salinity using gsw.O2sol_SP_pt with salinity and potential temperature ' \
        f'from dataset. Original oxygen concentration assumed to have been calculated using salinity = {ref_sal} PSU'
    return data

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
            _log.debug(a)
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
