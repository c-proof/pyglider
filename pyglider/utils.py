"""
Utilities that are used for processing scripts.
"""

import collections
import logging
from pathlib import Path

import gsw
import numpy as np
import xarray as xr
import yaml
from scipy.signal import argrelextrema
from scipy import signal
from datetime import datetime, date

from pyglider._version import __version__

_log = logging.getLogger(__name__)


def get_distance_over_ground(ds):
    """
    Add a distance over ground variable to a netcdf structure

    Parameters
    ----------
    ds : `xarray.Dataset`
        Must have variable ``latitude`` and ``longitude`` indexed
        by ``time`` dimension.

    Returns
    -------
    ds : `.xarray.Dataset`
        With ``distance_over_ground`` key.
    """

    good = ~np.isnan(ds.latitude + ds.longitude)
    if np.any(good):
        dist = gsw.distance(ds.longitude[good].values, ds.latitude[good].values) / 1000
        dist = np.roll(np.append(dist, 0), 1)
        dist = np.cumsum(dist)
        dist = np.interp(ds.time, ds.time[good], dist)
    else:
        dist = 0 * ds.latitude.values
    attr = {
        'long_name': 'distance over ground flown since mission start',
        'method': 'get_distance_over_ground',
        'units': 'km',
        'sources': 'latitude longitude',
    }
    ds['distance_over_ground'] = (('time'), dist, attr)
    return ds


def get_glider_depth(ds):
    """
    Get glider depth from pressure sensor.

    Parameters
    ----------
    ds : `xarray.Dataset`
        Must have variables ``pressure`` and ``latitude`` indexed
        by ``time`` dimension.  Assume pressure sensor in dbar.

    Returns
    -------
    ds : `.xarray.Dataset`
        With ``depth`` key.

    """
    good = np.where(~np.isnan(ds.pressure))[0]
    ds['depth'] = ds.pressure
    try:
        meanlat = ds.latitude.mean(skipna=True)
        ds['depth'].values = -gsw.z_from_p(
            ds.pressure.values, ds.latitude.fillna(meanlat).values
        )
    except AttributeError:
        pass
    # now we really want to know where it is, so interpolate:
    if len(good) > 0:
        ds['depth'].values = np.interp(
            np.arange(len(ds.depth)), good, ds['depth'].values[good]
        )

    attr = {
        'source': 'pressure',
        'long_name': 'glider depth',
        'standard_name': 'depth',
        'units': 'm',
        'comment': 'from science pressure and interpolated',
        'instrument': 'instrument_ctd',
        'observation_type': 'calulated',
        'accuracy': 1.0,
        'precision': 2.0,
        'resolution': 0.02,
        'platform': 'platform',
        'valid_min': 0.0,
        'valid_max': 2000.0,
        'reference_datum': 'surface',
        'positive': 'down',
    }
    ds['depth'].attrs = attr
    return ds


def get_profiles(ds, min_dp=10.0, inversion=3.0, filt_length=7, min_nsamples=14):
    """
    Not currently used...

    make two variables: profile_direction and profile_index; this version
    is good for lots of data.  Less good for sparse data
    """
    if 'pressure' not in ds:
        _log.warning(
            'No "pressure" variable in the data set; not searching for profiles'
        )
        return ds
    profile = ds.pressure.values * np.nan
    direction = ds.pressure.values * np.nan
    pronum = 1
    lastpronum = 0

    good = np.where(~np.isnan(ds.pressure))[0]
    p = np.convolve(
        ds.pressure.values[good], np.ones(filt_length) / filt_length, 'same'
    )
    dpall = np.diff(p)
    inflect = np.where(dpall[:-1] * dpall[1:] < 0)[0]
    for n, i in enumerate(inflect[:-1]):
        nprofile = inflect[n + 1] - inflect[n]
        inds = np.arange(good[inflect[n]], good[inflect[n + 1]] + 1) + 1
        dp = np.diff(ds.pressure[inds[[-1, 0]]])
        if (nprofile >= min_nsamples) and (np.abs(dp) > 10):
            _log.debug('Good')
            direction[inds] = np.sign(dp)
            profile[inds] = pronum
            lastpronum = pronum
            pronum += 1
        else:
            profile[good[inflect[n]] : good[inflect[n + 1]]] = lastpronum + 0.5

    attrs = collections.OrderedDict(
        [
            ('long_name', 'profile index'),
            ('units', '1'),
            ('comment', 'N = inside profile N, N + 0.5 = between profiles N and N + 1'),
            ('sources', 'time pressure'),
            ('method', 'get_profiles'),
            ('min_dp', min_dp),
            ('filt_length', filt_length),
            ('min_nsamples', min_nsamples),
        ]
    )
    ds['profile_index'] = (('time'), profile, attrs)

    attrs = collections.OrderedDict(
        [
            ('long_name', 'glider vertical speed direction'),
            ('units', '1'),
            ('comment', '-1 = ascending, 0 = inflecting or stalled, 1 = descending'),
            ('sources', 'time pressure'),
            ('method', 'get_profiles'),
        ]
    )
    ds['profile_direction'] = (('time'), direction, attrs)
    return ds


def get_profiles_new(ds, min_dp=10.0, filt_time=100, profile_min_time=300):
    """
    Find profiles in a glider timeseries:

    Parameters
    ----------
    ds : `xarray.Dataset`
        Must have *time* coordinate and *pressure* as a variable
    min_dp : float, default=10.0
        Minimum distance a profile must transit to be considered a profile, in dbar.
    filt_time : float, default=100
        Approximate length of time filter, in seconds.  Note that the filter
        is really implemented by sample, so the number of samples is
        ``filt_time / dt``
        where *dt* is the median time between samples in the time series.
    profile_min_time : float, default=300
        Minimum time length of profile in s.
    """

    if 'pressure' not in ds:
        _log.warning(
            'No "pressure" variable in the data set; not searching for profiles'
        )
        return ds

    profile = ds.pressure.values * 0
    direction = ds.pressure.values * 0
    pronum = 1

    good = np.where(np.isfinite(ds.pressure))[0]
    dt = float(
        np.median(np.diff(ds.time.values[good[:200000]]).astype(np.float64)) * 1e-9
    )
    _log.info(f'dt, {dt}')
    filt_length = int(filt_time / dt)

    min_nsamples = int(profile_min_time / dt)
    _log.info('Filt Len  %d, dt %f, min_n %d', filt_length, dt, min_nsamples)
    if filt_length > 1:
        p = np.convolve(
            ds.pressure.values[good], np.ones(filt_length) / filt_length, 'same'
        )
    else:
        p = ds.pressure.values[good]
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
    nmax = 0
    while (nmin < len(mins)) and (nmax < len(maxs)):
        nmax = np.where(maxs > mins[nmin])[0]
        if len(nmax) >= 1:
            nmax = nmax[0]
        else:
            break
        _log.debug(nmax)
        ins = range(int(mins[nmin]), int(maxs[nmax] + 1))
        _log.debug(f'{pronum}, {ins}, {len(p)}, {mins[nmin]}, {maxs[nmax]}')
        _log.debug(f'Down, {ins}, {p[ins[0]].values},{p[ins[-1]].values}')
        if (len(ins) > min_nsamples) and (
            np.nanmax(p[ins]) - np.nanmin(p[ins]) > min_dp
        ):
            profile[ins] = pronum
            direction[ins] = +1
            pronum += 1
        nmin = np.where(mins > maxs[nmax])[0]
        if len(nmin) >= 1:
            nmin = nmin[0]
        else:
            break
        ins = range(maxs[nmax], mins[nmin])
        _log.debug(f'{pronum}, {ins}, {len(p)}, {mins[nmin]}, {maxs[nmax]}')
        _log.debug(f'Up, {ins}, {p[ins[0]].values}, {p[ins[-1]].values}')
        if (len(ins) > min_nsamples) and (
            np.nanmax(p[ins]) - np.nanmin(p[ins]) > min_dp
        ):
            # up
            profile[ins] = pronum
            direction[ins] = -1
            pronum += 1

    attrs = collections.OrderedDict(
        [
            ('long_name', 'profile index'),
            ('units', '1'),
            ('comment', 'N = inside profile N, N + 0.5 = between profiles N and N + 1'),
            ('sources', 'time pressure'),
            ('method', 'get_profiles_new'),
            ('min_dp', min_dp),
            ('filt_length', filt_length),
            ('min_nsamples', min_nsamples),
        ]
    )
    ds['profile_index'] = (('time'), profile, attrs)

    attrs = collections.OrderedDict(
        [
            ('long_name', 'glider vertical speed direction'),
            ('units', '1'),
            ('comment', '-1 = ascending, 0 = inflecting or stalled, 1 = descending'),
            ('sources', 'time pressure'),
            ('method', 'get_profiles_new'),
        ]
    )
    ds['profile_direction'] = (('time'), direction, attrs)
    return ds


def get_derived_eos_raw(ds):
    """
    Calculate salinity, potential density, density, and potential temperature

    Parameters
    ----------
    ds : `xarray.Dataset`
        Must have *time* coordinate and *temperature*, *conductivity*, *pressure*,
        and *latitude* and *longitude* as variables.

    Returns
    -------
    ds : `xarray.Dataset`
        with *salinity*, *potential_density*, *density*, and *potential_temperature*
        as new variables.

    Notes
    -----
    Thermodynamic variables derived from the Gibbs seawater toolbox ``import gsw``.

    - *salinity*::
        gsw.conversions.SP_from_C(r, ds.temperature, ds.pressure)
    - *potential_density*::
        sa = gsw.SA_from_SP(ds['salinity'], ds['pressure'],
                            ds['longitude'], ds['latitude'])
        ct = gsw.CT_from_t(sa, ds['temperature'], ds['pressure'])
        ds['potential_density'] = (('time'),
                1000 + gsw.density.sigma0(sa, ct).values)
    - *density*::
        ds['density'] = (('time'), gsw.density.rho(
            ds.salinity, ds.temperature, ds.pressure).values)
    - *potential_temperature*::
        ds['potential_temperature'] = (('time'), gsw.conversions.pt0_from_t(
            ds.salinity, ds.temperature, ds.pressure).values)

    """

    # GPCTD and slocum ctd require a scale factor of 10 for conductivity.
    # Legato does not
    if 'S m' in ds.conductivity.units:
        r = 10 * ds.conductivity
    elif 'mS cm' in ds.conductivity.units:
        r = ds.conductivity
    else:
        raise ValueError(
            'Could not parse conductivity units in yaml. '
            "Expected 'S m-1' or 'mS cm-1'. "
            'Check yaml entry netcdf_variables: conductivity: units'
        )
    ds['salinity'] = (
        ('time'),
        gsw.conversions.SP_from_C(r, ds.temperature, ds.pressure).values,
    )
    attrs = collections.OrderedDict(
        [
            ('long_name', 'water salinity'),
            ('standard_name', 'sea_water_practical_salinity'),
            ('units', '1e-3'),
            ('comment', 'raw, uncorrected salinity'),
            ('sources', 'conductivity temperature pressure'),
            ('method', 'get_derived_eos_raw'),
            ('observation_type', 'calulated'),
            ('instrument', 'instrument_ctd'),
            ('valid_max', 40.0),
            ('valid_min', 0.0),
            ('accuracy', 0.01),
            ('precision', 0.01),
            ('resolution', 0.001),
        ]
    )
    attrs = fill_required_attrs(attrs)
    ds['salinity'].attrs = attrs
    long = ds.longitude.fillna(ds.longitude.mean(skipna=True))
    lat = ds.latitude.fillna(ds.latitude.mean(skipna=True))
    sa = gsw.SA_from_SP(ds['salinity'], ds['pressure'], long, lat)
    ct = gsw.CT_from_t(sa, ds['temperature'], ds['pressure'])
    ds['potential_density'] = (('time'), 1000 + gsw.density.sigma0(sa, ct).values)
    attrs = collections.OrderedDict(
        [
            ('long_name', 'water potential density'),
            ('standard_name', 'sea_water_potential_density'),
            ('units', 'kg m-3'),
            ('comment', 'raw, uncorrected salinity'),
            ('sources', 'salinity temperature pressure'),
            ('method', 'get_derived_eos_raw'),
            ('observation_type', 'calulated'),
            ('instrument', 'instrument_ctd'),
            ('accuracy', 0.01),
            ('precision', 0.01),
            ('resolution', 0.001),
        ]
    )
    attrs = fill_required_attrs(attrs)
    ds['potential_density'].attrs = attrs

    ds['density'] = (
        ('time'),
        gsw.density.rho(ds.salinity, ds.temperature, ds.pressure).values,
    )
    attrs = collections.OrderedDict(
        [
            ('long_name', 'Density'),
            ('standard_name', 'sea_water_density'),
            ('units', 'kg m-3'),
            ('comment', 'raw, uncorrected salinity'),
            ('observation_type', 'calulated'),
            ('sources', 'salinity temperature pressure'),
            ('instrument', 'instrument_ctd'),
            ('method', 'get_derived_eos_raw'),
            ('valid_min', 990.0),
            ('valid_max', 1040.0),
            ('accuracy', 0.01),
            ('precision', 0.01),
            ('resolution', 0.001),
        ]
    )
    attrs = fill_required_attrs(attrs)
    ds['density'].attrs = attrs
    ds['potential_temperature'] = (
        ('time'),
        gsw.conversions.pt0_from_t(ds.salinity, ds.temperature, ds.pressure).values,
    )
    attrs = collections.OrderedDict(
        [
            ('long_name', 'water potential temperature'),
            ('standard_name', 'sea_water_potential_temperature'),
            ('units', 'Celsius'),
            ('comment', 'raw, uncorrected salinity'),
            ('sources', 'salinity temperature pressure'),
            ('observation_type', 'calulated'),
            ('method', 'get_derived_eos_raw'),
            ('instrument', 'instrument_ctd'),
            ('accuracy', 0.002),
            ('precision', 0.001),
            ('resolution', 0.0001),
        ]
    )
    attrs = fill_required_attrs(attrs)
    ds['potential_temperature'].attrs = attrs

    return ds


def _time_to_datetime64(time):
    """
    Pass in a glider undecoded time (seconds since 1970-01-01), and
    get a np.datetime64[s] object back.
    """
    return time.astype('timedelta64[s]') + np.datetime64('1970-01-01')


def fill_required_attrs(attrs):
    required = {
        'comment': ' ',
        'accuracy': ' ',
        'precision': ' ',
        'platform': 'platform',
        'resolution': ' ',
        'ancillary_variables': ' ',
    }
    for k in required.keys():
        if k not in attrs.keys():
            attrs[k] = required[k]
    return attrs


def fill_required_qcattrs(attrs, varname):
    required = {
        'units': '1',
        'flag_values': np.array([1, 2, 3, 4, 9], dtype=np.int8),
        'valid_min': np.int8(1),
        'valid_max': np.int8(9),
        'flag_meanings': 'PASS NOT_EVALUATED SUSPECT FAIL MISSING',
        'standard_name': 'quality_flag',
        'long_name': 'Initial flag for {varname}',
    }
    for k in required.keys():
        if k not in attrs.keys():
            attrs[k] = required[k]
    return attrs


def get_file_id(ds):
    """
    Make a file id for a Dataset

    Parameters
    ----------
    ds : `xarray.Dataset`
        Dataset to make an id for.  The attributes of the Dataset
        must have *glider_name* and *glider_serial* in them.

    Return
    ------
    id : string
        Id = *glider_name* + *glider_serial* + "YYYYMMDDTHHMM"

    """

    _log.debug(ds.time)
    if not ds.time.dtype == 'datetime64[ns]':
        dt = _time_to_datetime64(ds.time.values[0])
    else:
        dt = ds.time.values[0].astype('datetime64[s]')
    _log.debug(f'dt, {dt}')
    id = (
        ds.attrs['glider_name']
        + ds.attrs['glider_serial']
        + '-'
        + dt.item().strftime('%Y%m%dT%H%M')
    )
    return id


def fill_metadata(ds, metadata, sensor_data):
    """
    Add metadata to a Dataset

    Parameters
    ----------
    ds : `xarray.Dataset`
        Dataset must have *longtidue*, *latitude*, and *time* values
    metadata : dict
        dictionary of attributes to add to the global attributes.  Usually
        taken from *deployment.yml* file.
    sensor_data : dict
        dictionary of device data to add to the global attributes.

    Returns
    -------
    ds : `xarray.Dataset`
        Dataset with attributes filled out.


    """
    good = ~np.isnan(ds.latitude.values + ds.longitude.values)
    if np.any(good):
        ds.attrs['geospatial_lat_max'] = np.max(ds.latitude.values[good])
        ds.attrs['geospatial_lat_min'] = np.min(ds.latitude.values[good])
        ds.attrs['geospatial_lon_max'] = np.max(ds.longitude.values[good])
        ds.attrs['geospatial_lon_min'] = np.min(ds.longitude.values[good])
    else:
        ds.attrs['geospatial_lat_max'] = np.nan
        ds.attrs['geospatial_lat_min'] = np.nan
        ds.attrs['geospatial_lon_max'] = np.nan
        ds.attrs['geospatial_lon_min'] = np.nan

    ds.attrs['geospatial_lat_units'] = 'degrees_north'
    ds.attrs['geospatial_lon_units'] = 'degrees_east'
    ds.attrs['netcdf_version'] = '4.0'  # TODO get this somehow...
    ds.attrs['history'] = (
        f'CPROOF glider toolbox pyglider version: {__version__}'
    )
    for k, v in metadata.items():
        ds.attrs[k] = v
    ds.attrs['featureType'] = 'trajectory'
    ds.attrs['cdm_data_type'] = 'Trajectory'
    ds.attrs['Conventions'] = 'CF-1.8'
    ds.attrs['standard_name_vocabulary'] = 'CF Standard Name Table v72'
    ds.attrs['date_created'] = str(np.datetime64('now')) + 'Z'
    ds.attrs['date_issued'] = str(np.datetime64('now')) + 'Z'
    ds.attrs['date_modified'] = ' '
    ds.attrs['id'] = get_file_id(ds)
    ds.attrs['title'] = ds.attrs['id']

    dt = ds.time.values
    ds.attrs['time_coverage_start'] = '%s' % dt[0]
    ds.attrs['time_coverage_end'] = '%s' % dt[-1]
    
    # make sure this is ISO readable....
    ds.attrs['deployment_start'] = str(dt[0])[:19]
    ds.attrs['deployment_end'] = str(dt[-1])[:19]

    ds.attrs['processing_level'] = (
        'Level 0 (L0) processed data timeseries; no corrections or data screening'
    )

    for k, v in sensor_data.items():
        ds.attrs[k] = str(v)
    ds.attrs = collections.OrderedDict(sorted(ds.attrs.items()))

    return ds


def _zero_screen(val):
    val[val == 0] = np.nan
    return val


def nmea2deg(nmea):
    """
    Convert a NMEA float to a decimal degree float.  e.g. -12640.3232 = -126.6721
    """
    deg = np.fix(nmea / 100) + np.sign(nmea) * np.remainder(np.abs(nmea), 100) / 60
    return deg


def oxygen_concentration_correction(ds, ncvar):
    """
    Correct oxygen signal for salinity signal

    Parameters
    ----------
    ds : `xarray.Dataset`
        Should have *oxygen_concentration*, *potential_temperature*, *salinity*,
        on a *time* coordinate.
    ncvar : dict
        dictionary with netcdf variable definitions in it.  Should have
        *oxygen_concentration* as a key, which itself should specify
        a *reference_salinity* and have *correct_oxygen* set to ``"True"``.

    Returns
    -------
    ds : `xarray.Dataset`
        With *oxygen_concentration* corrected for the salinity effect.
    """

    oxy_yaml = ncvar['oxygen_concentration']
    if 'reference_salinity' not in oxy_yaml.keys():
        _log.warning(
            'No reference_salinity found in oxygen deployment yaml. '
            'Assuming reference salinity of 0 psu'
        )
        ref_sal = 0
    else:
        ref_sal = float(oxy_yaml['reference_salinity'])
    _log.info(f'Correcting oxygen using reference salinity {ref_sal} PSU')
    ds_oxy = ds.oxygen_concentration[~np.isnan(ds.oxygen_concentration)]
    # Match the nearest temperature and salinity values from their timestamps
    ds_temp = ds.potential_temperature[~np.isnan(ds.potential_temperature)].reindex(
        time=ds_oxy.time, method='nearest'
    )
    ds_sal = ds.salinity[~np.isnan(ds.salinity)].reindex(
        time=ds_oxy.time, method='nearest'
    )
    o2_sol = gsw.O2sol_SP_pt(ds_sal, ds_temp)
    o2_sat = ds_oxy / gsw.O2sol_SP_pt(ds_sal * 0 + ref_sal, ds_temp)
    ds['oxygen_concentration'].values[~np.isnan(ds.oxygen_concentration)] = (
        o2_sat * o2_sol
    )
    ds['oxygen_concentration'].attrs['oxygen_concentration_QC:RTQC_methodology'] = (
        f'oxygen concentration corrected for salinity using gsw.O2sol_SP_pt with'
        f'salinity and potential temperature from dataset. Original oxygen '
        f'concentration assumed to have been calculated using salinity = '
        f'{ref_sal} PSU'
    )
    return ds


def bar2dbar(val):
    """
    convert val bar to dbar
    """
    return val * 10


def dbar2bar(val):
    """
    convert val dbar to bar
    """
    return val / 10


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
        if (
            len(ind) > 0
            and len(ind) < (ind[-1] - ind[0])
            and len(ind) > (ind[-1] - ind[0]) * 0.05
        ):
            int = np.arange(ind[0], ind[-1])
            data[:, j][ind[0] : ind[-1]] = np.interp(int, ind, data[ind, j])
    return data


def find_gaps(sample_time, timebase, maxgap):
    """
    Return an index into *timebase* where True are times in gaps of *sample_time* larger
    than maxgap.
    """
    # make sure sample times are strictly increasing
    sample_time = np.sort(sample_time)

    # figure out which sample each time in time base belongs to:
    time_index = np.searchsorted(sample_time, timebase, side='right')
    time_index = np.clip(time_index, 0, len(sample_time) - 1)

    # figure out the space between sample pairs
    dt = np.concatenate(([0], np.diff(sample_time)))
    # get the gap size for each timebase data point:
    ddt = dt[time_index]

    # get the indices of timebase that are too large and account for the
    # degenerate case when a timebase point falls directly on a sample time.
    index = ~np.logical_or((ddt <= maxgap), (np.isin(timebase, sample_time)))

    return index


def _parse_gliderxml_pos(fname):
    """
    DEPRECATED: use slocum.parse_gliderState instead

    returns lon, lat, timestr
    """

    xmln = fname
    from bs4 import BeautifulSoup

    with open(xmln, 'r') as fin:
        y = BeautifulSoup(fin, features='xml')
        time = None
        lon = None
        lat = None
        for a in y.find_all('valid_location'):
            try:
                dtime = np.datetime64(a.time.text)
                if dtime > np.datetime64('2019-07-19'):
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


def _parse_gliderxml_surfacedatetime(fname):
    """
    DEPRECATED: use slocum.parse_gliderState instead

    returns lon, lat, timestr
    """

    xmln = fname
    from bs4 import BeautifulSoup

    with open(xmln, 'r') as fin:
        y = BeautifulSoup(fin, features='xml')
        time = None
        for a in y.find_all('report'):
            _log.debug(a)
            if a.text is not None:
                try:
                    time = np.append(time, np.datetime64(a.text))
                except:
                    pass
    return time


def example_gridplot(
    filename,
    outname,
    toplot=['potential_temperature', 'salinity', 'oxygen_concentration'],
    pdenlevels=np.arange(10, 30, 0.5),
    dpi=200,
    ylim=None,
):
    """
    Smoketest for plotting the netcdf files.
    """

    import matplotlib.pyplot as plt

    ntoplot = len(toplot)
    with xr.open_dataset(filename) as ds:
        fig, axs = plt.subplots(
            nrows=ntoplot,
            constrained_layout=True,
            figsize=(7, 3 * ntoplot),
            sharex=True,
            sharey=True,
        )
        for ax, vname in zip(axs, toplot):
            ds[vname].plot.pcolormesh(ax=ax)
            (ds['potential_density'] - 1000).plot.contour(ax=ax, levels=pdenlevels)
            if ylim:
                ax.set_ylim(ylim)
        fig.savefig(outname, dpi=dpi)


def _get_deployment(deploymentyaml):
    """
    Take the list of files in *deploymentyaml* and parse them
    for deployment information, with subsequent files overwriting
    previous files.
    """
    if isinstance(deploymentyaml, str):
        deploymentyaml = [
            deploymentyaml,
        ]
    deployment = {}
    for nn, d in enumerate(deploymentyaml):
        with open(d) as fin:
            deployment_ = yaml.safe_load(fin)
            for k in deployment_:
                deployment[k] = deployment_[k]

    return deployment


def _any_newer(dirname, filename):
    """
    Check if any files in dirname are newer than filename
    """
    filename = Path(filename)
    dirname = Path(dirname)
    print(filename, filename.exists())
    if not filename.exists():
        return True

    mod_time = filename.stat().st_mtime
    is_newer = False
    for file_path in dirname.iterdir():
        if file_path.is_file():
            if file_path.stat().st_mtime > mod_time:
                is_newer = True
                break

    return is_newer


def _get_glider_name_slocum(current_directory):
    glider = current_directory.parts[-2]
    mission = current_directory.parts[-1]
    print(f'Glider {glider} and mission: {mission}')
    slocum_glider = glider[4:]
    if slocum_glider[-4:-3].isnumeric():
        slocum_glider = slocum_glider[:-4] + '_' + slocum_glider[-4:]
    else:
        slocum_glider = slocum_glider[:-3] + '_' + slocum_glider[-3:]

    if slocum_glider == 'walle_652':
        slocum_glider = 'wall_e_652'
    return glider, mission, slocum_glider

def flag_conductivity_in_depth_space(ts0, d_profile=50, dz=5, clean_stdev=3, accuracy=None):
    """
    Flag conductivity as QC1 (good) or QC4 (bad) using profile bins and depth bins.

    Conductivity data are grouped into profile bins of width `d_profile` and depth bins
    of width `dz`. Within each depth bin, points farther than 5 standard deviations from
    the mean are temporarily excluded. A second mean and standard deviation are then
    computed from the remaining points, and values farther than `clean_stdev` standard
    deviations from that cleaned mean are flagged as QC4. If `accuracy` is provided,
    deviations smaller than `accuracy` are not flagged.

    Parameters
    ----------
    ts0 : xarray.Dataset
        Timeseries dataset containing conductivity, depth, and profile_index.

    d_profile : float, optional
        Width of the profile bins.

    dz : float, optional
        Width of the depth bins in meters.

    clean_stdev : float, optional
        Number of standard deviations used in the second-pass flagging step.

    accuracy : float or None, optional
        Sensor accuracy threshold. Deviations smaller than this are not flagged.

    Returns
    -------
    qc : np.ndarray
        Array of QC flags with the same shape as ts0.conductivity.
        Good data are flagged as 1, bad data as 4.
    """
    ts = ts0.copy(deep=True).load()
    ts = ts.where(np.isfinite(ts.conductivity), drop=False)
    _log.info('Flagging conductivity in profile-depth space with d_profile=%s, dz=%s, clean_stdev=%s, accuracy=%s', d_profile, dz, clean_stdev, accuracy)

    prof_vals = ts.profile_index.values
    depth_vals = ts.depth.values
    cond_vals = ts.conductivity.values

    prof_bins = np.arange(np.nanmin(prof_vals), np.nanmax(prof_vals) + d_profile, d_profile)
    zbins = np.arange(np.nanmin(depth_vals), np.nanmax(depth_vals) + dz, dz)

    qc = np.full(cond_vals.shape, 4, dtype=int)
    finite_cond = np.isfinite(cond_vals)

    acc_thresh = 0 if accuracy is None else accuracy

    for n in range(len(prof_bins) - 1):
        ind_profbin = (
            (prof_vals >= prof_bins[n]) &
            (prof_vals < prof_bins[n + 1]) &
            finite_cond
        )

        if not np.any(ind_profbin):
            continue

        cond = cond_vals[ind_profbin]
        depth = depth_vals[ind_profbin]

        ind_bad_z = np.zeros(len(cond), dtype=bool)

        for m in range(len(zbins) - 1):
            ind_zbin = (depth >= zbins[m]) & (depth < zbins[m + 1])

            if not np.any(ind_zbin):
                continue

            var_z = cond[ind_zbin]
            var_mean = np.nanmean(var_z)
            var_std = np.nanstd(var_z)

            if not np.isfinite(var_std) or var_std == 0:
                ind_bad = np.zeros_like(var_z, dtype=bool)
            else:
                ind_flag = (
                    (np.abs(var_z - var_mean) > 5 * var_std) &
                    (np.abs(var_z - var_mean) > acc_thresh)
                )

                if np.all(ind_flag):
                    ind_bad = ind_flag
                else:
                    clean_mean = np.nanmean(var_z[~ind_flag])
                    clean_std = np.nanstd(var_z[~ind_flag])

                    if not np.isfinite(clean_std) or clean_std == 0:
                        ind_bad = np.zeros_like(var_z, dtype=bool)
                    else:
                        ind_bad = (
                            (np.abs(var_z - clean_mean) > clean_stdev * clean_std) &
                            (np.abs(var_z - clean_mean) > acc_thresh)
                        )

            ind_bad_z[ind_zbin] = ind_bad

        qc_subset = np.where(ind_bad_z, 4, 1)
        qc[ind_profbin] = qc_subset

    return qc

def interpolate_over_salinity_NANs(ds): 
    """
    Function applied to the dataset before finding the internal temperature. Function interpolates temperature over bad data and small data gaps
    to prevent errors from affecting the neighbiuring cells. 

    Parameters
    ----------
    ds: DataArray
        Timeseries of mission data 
    
    Returns
    ----------
    interp: DataArray
        Timeseries of interpolated temperature
    
    """
    _log.info('Interpolating temperature over salinity NaNs and small data gaps before applying thermal lag correction')
    interp = ds["temperature"].where(ds["temperature_QC"] != 4)
    qc4 = (ds["temperature_QC"] == 4)
    qc4_buf = qc4.rolling(time=5, center=True, min_periods=1).max().astype(bool)
    interp = interp.where(~qc4_buf)

    interp = interp.interpolate_na(
        dim="time",
        method="linear",
        max_gap=np.timedelta64(60, "s"))

    return interp

def apply_thermal_lag(ds, fn, alpha, tau, interpolate_filter = None):
    """
    Function from Garau et al. (2011): estimates temperature inside the conductivity cell then recaluclates salinity

    Parameters
    ----------
    ds: DataArray
        Timeseries of mission data 

    fn: float 
        Sampling frequency of the sensor

    alpha : float
        Thermal lag strength constant for the sensor. 

    tau: float
        Thermal lag time constant for the sensor. 

    interpolate_filter: callable or None, optional
        Function applied to the dataset before finding the internal temperature. Function interpolates over bad data and small data gaps
        to prevent errors from affecting the neighbiuring cells. 
    
    Returns 
    ----------
    sal: DataArray
        Timeseries of salinity_adjusted calculated using the internal temperature of the conductivity cell.
    """
    if interpolate_filter is not None: 
        temp = interpolate_filter(ds)
        _log.info('Interpolating over bad data and small data gaps before applying thermal lag correction')

    else: 
        temp = ds.temperature

    _log.info('Applying thermal lag correction with alpha = %s, tau = %s, and sampling frequency = %s', alpha, tau, fn)

    a = 4 * fn * alpha * tau / (1 + 4*fn*tau)
    b = 1 - 2 * a / alpha
    aa = [1, b]
    bb = [a, -a]
    tempcorr = temp.values.copy()
    tempcell = temp.values.copy()
    good = ~np.isnan(tempcell)
    tempcorr[good] = signal.lfilter(bb, aa, temp.values[good])
    tempcell = tempcell - tempcorr
    sal = gsw.SP_from_C(ds.conductivity * 10, tempcell, ds.pressure)

    return sal

def flag_CTD_data(ts0, clean_stdev=3, accuracy=None):
    """
    Wrapper function to flag CTD data.

    Uses `flag_conductivity_in_depth_space` to flag conductivity as QC1 (good)
    or QC4 (bad) in profile-depth space. Conductivity and salinity are then
    flagged as QC4 wherever conductivity is flagged as QC4.

    Creates `conductivity_QC`, `salinity_QC`, and `temperature_QC` if they do
    not already exist.

    Parameters
    ----------
    ts0 : xarray.Dataset
        Timeseries of mission data.

    clean_stdev : float, optional
        Number of standard deviations from the cleaned mean for data to be
        flagged as QC4.

    Returns
    -------
    ts : xarray.Dataset
        Timeseries of mission data with `conductivity_QC`, `salinity_QC`,
        and `temperature_QC`.
    """
    _log.info('Sceening CTD data')

    ts = ts0.copy()

    ts["conductivity"] = ts["conductivity"].where(ts["conductivity"] >= 0.1)

    cond_qc = flag_conductivity_in_depth_space(
        ts,
        d_profile=50,
        dz=5,
        clean_stdev=clean_stdev,
        accuracy=accuracy
    )

    if "conductivity_QC" not in ts.data_vars:
        _log.info('Adding conductivity_QC variable to dataset')

        ts["conductivity_QC"] = xr.DataArray(
            np.ones(ts["conductivity"].shape, dtype=int),
            dims=ts["conductivity"].dims,
            coords=ts["conductivity"].coords,
        )

    if "salinity_QC" not in ts.data_vars:
        _log.info('Adding salinity_QC variable to dataset')
        ts["salinity_QC"] = xr.DataArray(
            np.ones(ts["salinity"].shape, dtype=int),
            dims=ts["salinity"].dims,
            coords=ts["salinity"].coords,
        )

    if "temperature_QC" not in ts.data_vars:
        _log.info('Adding temperature_QC variable to dataset')
        ts["temperature_QC"] = xr.DataArray(
            np.ones(ts["temperature"].shape, dtype=int),
            dims=ts["temperature"].dims,
            coords=ts["temperature"].coords,
        )

    ts["conductivity_QC"] = xr.where(cond_qc == 4, 4, ts["conductivity_QC"])
    ts["salinity_QC"] = xr.where(ts["conductivity_QC"] == 4, 4, ts["salinity_QC"])

    return ts

def adjust_CTD(ts, deploymentyaml, alpha=None, tau=None, dTdC=None, interpolate_filter=None):
    """
    Pulls correction constants from `deploymentyaml`. If `alpha`, `tau`, or `dTdC`
    differ from the values in the YAML file, the values provided as function arguments
    are used and a warning is issued.
    
    Applies conductivity–temperature lag correction and thermal lag correction when
    the corresponding constants are not `None` or 0. This produces the variables
    `temperature_adjusted` and `salinity_adjusted`. The variables
    `potential_density_adjusted` and `potential_temperature_adjusted` are derived
    from the adjusted temperature and salinity.
    
    Parameters
    ----------
    ts : xarray.Dataset
        Time series of mission data.
    
    deploymentyaml : str or list
        Path to a YAML file containing deployment information for the glider.
    
        If a list is provided, YAML files are read in order, and top-level keys
        in later files overwrite those in earlier files.
    
    alpha : float, optional
        Thermal lag correction parameter alpha. Default is None.
    
    tau : float, optional
        Thermal lag correction parameter tau. Default is None.
    
    dTdC : float, optional
        Time lag (seconds) between temperature and conductivity sensors. Default is None.
    
    interpolate_filter: callable or None, optional
        Function applied to the dataset before finding the internal temperature. Function interpolates over bad data and small data gaps
        to prevent errors from affecting the neighbiuring cells. Default is None.
    Returns
    -------
    ts : xarray.Dataset
        Time series dataset with the additional variables:
        `temperature_adjusted`, `salinity_adjusted`,
        `potential_density_adjusted`, and `potential_temperature_adjusted`.
        Metadata are updated to reflect applied corrections.
    """
    logger = logging.getLogger(__name__)
    _log.info('Adjusting CTD data')

    atr = deploymentyaml.get("glider_devices", {}).get("ctd", {})
    thermal = atr.get("Thermal_lag_constants_[alpha,tau]")
    _log.info('CTD thermal lag constants from YAML: %s', thermal)
    yaml_vals = {
        "alpha": thermal[0] if thermal and len(thermal) > 0 else None,
        "tau":   thermal[1] if thermal and len(thermal) > 1 else None,
        "dTdC":  atr.get("dTdC"),
    }

    kw_vals = {
        "alpha": alpha,
        "tau": tau,
        "dTdC": dTdC,
    }

    out = {}
    for key in yaml_vals:
        y = yaml_vals[key]
        k = kw_vals[key]

        if k is not None:
            if y is not None and y != k and logger is not None:
                logger.warning(
                    "%s differs between YAML (%r) and kwargs (%r); using kwargs.",
                    key, y, k
                )
            out[key] = k
        else:
            out[key] = y

    alpha = out.get("alpha", {})
    tau = out.get("tau", {})
    dTdC = out.get("dTdC", {})

    if all(out.get(k) is None for k in ["alpha", "tau", "dTdC"]):
        raise ValueError(
            "Missing required CTD constants after checking kwargs and YAML: alpha, tau, dTdC"
        )

    temp_adj = ts.temperature.copy()
    temp_adj.attrs = ts.temperature.attrs.copy()
    temp_adj.attrs["comment"] = "temperature [degC]"

    if dTdC not in (None, 0):
        _log.info('Interpolating temperature data forward by %s seconds', dTdC)
        dt = np.timedelta64(dTdC, "s")
        temp_adj = temp_adj.interp(time=ts.time + dt)

        temp_adj.attrs["history"] = "temperature [degC] adjusted by CT lag"
        temp_adj.attrs["time_lag"] = f"{dTdC} second CT lag corrected"
        ts.attrs["dTdC"] = f"{dTdC} second CT lag corrected"
    else:
        temp_adj.attrs["comment"] = "equivalent to raw temperature"
        ts.attrs["dTdC"] = "No CT lag applied"

    ts["temperature_adjusted"] = temp_adj

    if tau is not None:
        dt = np.diff(ts.time.values).astype("timedelta64[s]").astype(int)
        vals, counts = np.unique(dt, return_counts=True)
        srate = vals[np.argmax(counts)]

        fs = 1 / float(srate)
        fn = 0.5 * fs

        s = apply_thermal_lag(ts, fn, alpha=alpha, tau=tau, interpolate_filter=interpolate_filter)

        sal_adj = xr.where(ts.salinity_QC == 1, s, ts.salinity)
        sal_adj.attrs = ts.salinity.attrs.copy()
        sal_adj.attrs["history"] = (
            f"adjusted salinity [psu] using thermal lag correction "
            f"(alpha={alpha}, tau={tau})"
        )

        if dTdC not in (None, 0):
            sal_adj.attrs["sources"] = (
                f"conductivity pressure temperature_adjusted "
                f"(corrected for {dTdC} second CT lag)"
            )

        ts["salinity_adjusted"] = sal_adj
        ts.attrs['correction_constants_alpha'] = alpha 
        ts.attrs['correction_constants_tau'] = tau 

    else:
        _log.info('No thermal lag correction applied; calculating salinity_adjusted using temperature_adjusted and raw conductivity')
        sal_adj = xr.DataArray(
            gsw.conversions.SP_from_C(
                10 * ts["conductivity"],
                ts["temperature_adjusted"],

                ts.pressure,
            ).values,
            dims=ts.salinity.dims,
            coords=ts.salinity.coords,
        )
        ts.attrs['correction_constants_alpha'] = "None"
        ts.attrs['correction_constants_tau'] =  "None"

        sal_adj.attrs = ts.salinity.attrs.copy()

        if dTdC is not None:
            sal_adj.attrs["time_lag"] = "found using temperature_adjusted"

        ts["salinity_adjusted"] = sal_adj

    ts["salinity_adjusted"].attrs = sal_adj.attrs

    ts.attrs["quality_flags"] = (
        "1 = good data; 3 = bad data, potentially correctable; "
        "4 = bad data; 8 = estimated data"
    )

    ts["conductivity"].attrs["comment"] = "raw conductivity"
    ts["conductivity_QC"] = ts.conductivity_QC

    ts["temperature"].attrs["comment"] = "raw temperature [degC]"
    ts["temperature_QC"] = ts.temperature_QC

    ts["temperature_adjusted_QC"] = ts["temperature_QC"]

    ts["salinity"].attrs["comment"] = "raw salinity [psu]"
    ts["salinity_adjusted_QC"] = ts["salinity_QC"]

    ts["density"].attrs["comment"] = "raw density"
    ts["density_QC"] = ts["salinity_QC"]

    ts["potential_density"].attrs["history"] = (
        "calculated using raw salinity and temperature"
    )
    ts["potential_temperature"].attrs["history"] = (
        "calculated using raw salinity and temperature"
    )

    long = ts.longitude.fillna(ts.longitude.mean(skipna=True))
    lat = ts.latitude.fillna(ts.latitude.mean(skipna=True))

    sa_adj = gsw.SA_from_SP(ts["salinity_adjusted"], ts["pressure"], long, lat)
    ct_adj = gsw.CT_from_t(sa_adj, ts["temperature_adjusted"], ts["pressure"])

    _log.info('Calculating potential density and potential temperature using adjusted salinity and temperature')
    ts["potential_density_adjusted"] = (
        ("time"), 1000 + gsw.density.sigma0(sa_adj, ct_adj).values
    )

    ts["potential_density"].attrs = ts.potential_density.attrs.copy()

    ts["potential_density_adjusted"].attrs["history"] = (
        "calculated using salinity_adjusted and temperature_adjusted"
    )
    ts["potential_density_adjusted"].attrs["sources"] = (
        "salinity_adjusted temperature_adjusted pressure"
    )

    ts["potential_density_adjusted_QC"] = ts["salinity_adjusted_QC"]
    ts["potential_density_adjusted_QC"].attrs = ts.salinity_adjusted_QC.attrs.copy()

    ts["potential_temperature_adjusted"] = (
        ("time"),
        gsw.conversions.pt0_from_t(
            ts.salinity_adjusted,
            ts.temperature_adjusted,
            ts.pressure,
        ).values,
    )

    ts["potential_temperature_adjusted"].attrs = (
        ts.potential_temperature.attrs.copy()
    )
    ts["potential_temperature_adjusted"].attrs["history"] = (
        "calculated using salinity_adjusted and temperature_adjusted"
    )
    ts["potential_temperature_adjusted"].attrs["sources"] = (
        "salinity_adjusted temperature_adjusted pressure"
    )

    ts["potential_temperature_adjusted_QC"] = ts["salinity_adjusted_QC"]
    ts["potential_temperature_adjusted_QC"].attrs = (
        ts.salinity_adjusted_QC.attrs.copy()
    )

    processing_date = date.today().strftime("%Y%m%d")

    vars_ = [
        "salinity_adjusted",
        "temperature_adjusted",
        "potential_density_adjusted",
        "potential_temperature_adjusted",
    ]

    for var in vars_:
        ts[var].attrs["processing_date"] = processing_date

    QC_COMMENT = (
        "1 = good data; 3 = bad data, potentially correctable; "
        "4 = bad data; 8 = estimated data"
    )

    for k in ts.data_vars:
        if k.endswith("_QC"):
            ts[k].attrs["comment"] = QC_COMMENT

    return ts

def maskQC4(ds):
    """
    Optional: 
    Masks QC4 samples in data variables (set to NaN) so gridding ignores them. Only QC1 (good) data 
    are gridded
    
    Parameters
    ----------
    ds: DataArray
        Timeseries of a data   

    Returns
    ----------
    ds: DataArray
        Timeseries of a data with QC4 data masked    

    
    """
    _log = logging.getLogger(__name__)

    ds = ds.copy() 
    _log.info('Masking QC4 data in dataset')
    for k in list(ds.data_vars):
        # skip QC variables themselves
        if k.endswith("_QC"):
            continue

        qc_name = f"{k}_QC"
        if qc_name not in ds:
            continue

        # mask data where QC == 4, preserving dims/coords
        ds[k] = ds[k].where(ds[qc_name] != 4)
        ds[qc_name] = ds[qc_name].where(ds[qc_name] != 4)
        
    return ds

def interpolate_vertical(var, attr):
    """
    Optional: Interpolates variables over gaps of 50m. If the variable has 'QC_protocol' in the metadata, NaN gaps are filled
    with QC1 (good flag) 

    Parameters
    ----------
    var: DataArray
        Timeseries of a data variable

    attr: 
    
    """ 
    _log.info('Interpolating variable %s over depth with max gap of 50m', var.name)
    # QC variables: fill interpolatable NaN gaps with 1
    if 'QC_protocol' in attr.attrs.values():
        interp = var.interpolate_na(dim="depth", method="nearest", max_gap=50)
        filled = np.isnan(var) & np.isfinite(interp)
        return xr.where(filled, 1, var)

    # Continuous variables: linear interpolation
    return var.interpolate_na(dim="depth", method="linear", max_gap=50)

__all__ = [
    'get_distance_over_ground',
    'get_glider_depth',
    'get_profiles_new',
    'get_derived_eos_raw',
    'fill_metadata',
    'nmea2deg',
    'gappy_fill_vertical',
    'oxygen_concentration_correction',
    'flag_conductivity_in_depth_space',
    'interpolate_over_salinity_NANs'
    'apply_thermal_lag',
    'flag_CTD_data',
    'adjust_CTD',
    'maskQC4',
    'interpolate_vertical'
]