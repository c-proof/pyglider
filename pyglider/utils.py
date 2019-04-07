import collections
import seawater
import xarray as xr
import numpy as np


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
    ds['depth'] = seawater.eos80.dpth(ds.pressure, ds.latitude.mean())
    # now we really want to know where it is, so interpolate:
    if len(good) > 0:
        ds['depth'] = np.interp(np.arange(len(ds.depth)),
                                good, ds['depth'][good])

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
    make two variables: profile_direction and profile_index
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
        inds = np.arange(good[inflect[n]], good[inflect[n+1]]+1)
        dp = np.diff(ds.pressure[inds[[-1, 0]]])
        if ((nprofile >= min_nsamples) and (np.abs(dp) > 10)):

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


def get_derived_eos_raw(ds):
    r = 10 * ds.conductivity / seawater.constants.c3515
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
    Pass in a glider undecodeed time (seconds since 1970-01-01), and
    get a np.datetime64[s] object back.
    """
    return (time.astype('timedelta64[s]') +
            np.datetime64('1970-01-01'))


def fill_required_attrs(attrs):
    print(attrs)
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
    dt = time_to_datetime64(ds.time.values)
    id = (ds.attrs['glider_name'] + ds.attrs['glider_serial'] + '-' +
                      dt[0].item().strftime('%Y%m%dT%H%M'))
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
    ds.attrs['title'] = ds.attrs['id']

    dt = time_to_datetime64(ds.time.values)
    ds.attrs['time_coverage_start'] = '%s' % dt[0]
    ds.attrs['time_coverage_end'] = '%s' % dt[-1]

    ds.attrs['processing_level'] = 'Level 1 (L1) processed data timeseries; no corrections or data screening'

    ds.attrs = collections.OrderedDict(sorted(ds.attrs.items()))

    return ds

def _zero_screen(val):
    val[val==0] = np.NaN
    return val


def nmea2deg(nmea):
    return np.fix(nmea / 100) + np.remainder(nmea, 100) / 60


def bar2dbar(val):
    return val * 10.0


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
        if len(ind) > 0 and len(ind) < (ind[-1] - ind[0]):
            int = np.arange(ind[0], ind[-1])
            data[:, j][ind[0]:ind[-1]] = np.interp(int, ind, data[ind, j])
    return data
