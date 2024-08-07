"""
Routines that are used for common processing of netcdf files after they have
been converted to standard timeseries.
"""
import logging
import xarray as xr
import numpy as np
import pyglider.utils as utils
import os
import yaml
import netCDF4
import scipy.stats as stats

_log = logging.getLogger(__name__)


def extract_timeseries_profiles(inname, outdir, deploymentyaml, force=False):
    """
    Extract and save each profile from a timeseries netCDF.

    Parameters
    ----------
    inname : str or Path
        netcdf file to break into profiles

    outdir : str or Path
        directory to place profiles

    deploymentyaml : str or Path
        location of deployment yaml file for the netCDF file.  This should
        be the same yaml file that was used to make the timeseries file.

    force : bool, default False
        Force an overwite even if profile netcdf already exists
    """
    try:
        os.mkdir(outdir)
    except FileExistsError:
        pass

    deployment = utils._get_deployment(deploymentyaml)

    meta = deployment['metadata']
    with xr.open_dataset(inname) as ds:
        _log.info('Extracting profiles: opening %s', inname)
        profiles = np.unique(ds.profile_index)
        profiles = [p for p in profiles if (~np.isnan(p) and not (p % 1)
                                            and (p > 0))]
        for p in profiles:
            ind = np.where(ds.profile_index == p)[0]
            dss = ds.isel(time=ind)
            outname = outdir + '/' + utils.get_file_id(dss) + '.nc'
            _log.info('Checking %s', outname)
            if force or (not os.path.exists(outname)):
                # this is the id for the whole file, not just this profile..
                dss['trajectory'] = utils.get_file_id(ds).encode()
                trajlen = len(utils.get_file_id(ds).encode())
                dss['trajectory'].attrs['cf_role'] = 'trajectory_id'
                dss['trajectory'].attrs['comment'] = (
                    'A trajectory is a single'
                    'deployment of a glider and may span multiple data files.')
                dss['trajectory'].attrs['long_name'] = 'Trajectory/Deployment Name'

                # profile-averaged variables....
                profile_meta = deployment['profile_variables']
                if 'water_velocity_eastward' in dss.keys():
                    dss['u'] = dss.water_velocity_eastward.mean()
                    dss['u'].attrs = profile_meta['u']

                    dss['v'] = dss.water_velocity_northward.mean()
                    dss['v'].attrs = profile_meta['v']
                elif 'u' in profile_meta:
                    dss['u'] = profile_meta['u'].get('_FillValue', np.nan)
                    dss['u'].attrs = profile_meta['u']

                    dss['v'] = profile_meta['v'].get('_FillValue', np.nan)
                    dss['v'].attrs = profile_meta['v']
                else:
                    dss['u'] = np.nan
                    dss['v'] = np.nan


                dss['profile_id'] = np.int32(p)
                dss['profile_id'].attrs = profile_meta['profile_id']
                if '_FillValue' not in dss['profile_id'].attrs:
                    dss['profile_id'].attrs['_FillValue'] = -1
                dss['profile_id'].attrs['valid_min'] = np.int32(dss['profile_id'].attrs['valid_min'])
                dss['profile_id'].attrs['valid_max'] = np.int32(dss['profile_id'].attrs['valid_max'])

                dss['profile_time'] = dss.time.mean()
                dss['profile_time'].attrs = profile_meta['profile_time']
                # remove units so they can be encoded later:
                try:
                    del dss.profile_time.attrs['units']
                    del dss.profile_time.attrs['calendar']
                except KeyError:
                    pass
                dss['profile_lon'] = dss.longitude.mean()
                dss['profile_lon'].attrs = profile_meta['profile_lon']
                dss['profile_lat'] = dss.latitude.mean()
                dss['profile_lat'].attrs = profile_meta['profile_lat']

                dss['lat'] = dss['latitude']
                dss['lon'] = dss['longitude']
                dss['platform'] = np.int32(1)
                comment = (meta['glider_model'] + ' operated by ' +
                           meta['institution'])
                dss['platform'].attrs['comment'] = comment
                dss['platform'].attrs['id'] = (
                    meta['glider_name'] + meta['glider_serial'])
                dss['platform'].attrs['instrument'] = 'instrument_ctd'
                dss['platform'].attrs['long_name'] = (
                    meta['glider_model'] + dss['platform'].attrs['id'])
                dss['platform'].attrs['type'] = 'platform'
                dss['platform'].attrs['wmo_id'] = meta['wmo_id']
                if '_FillValue' not in dss['platform'].attrs:
                    dss['platform'].attrs['_FillValue'] = -1


                dss['lat_uv'] = np.nan
                dss['lat_uv'].attrs = profile_meta['lat_uv']
                dss['lon_uv'] = np.nan
                dss['lon_uv'].attrs = profile_meta['lon_uv']
                dss['time_uv'] = np.nan
                dss['time_uv'].attrs = profile_meta['time_uv']

                dss['instrument_ctd'] = np.int32(1.0)
                dss['instrument_ctd'].attrs = profile_meta['instrument_ctd']
                if '_FillValue' not in dss['instrument_ctd'].attrs:
                    dss['instrument_ctd'].attrs['_FillValue'] = -1

                dss.attrs['date_modified'] = str(np.datetime64('now')) + 'Z'

                # ancillary variables: link and create with values of 2.  If
                # we dont' want them all 2, then create these variables in the
                # time series
                to_fill = ['temperature', 'pressure', 'conductivity',
                        'salinity', 'density', 'lon', 'lat', 'depth']
                for name in to_fill:
                    qcname = name + '_qc'
                    dss[name].attrs['ancillary_variables'] = qcname
                    if qcname not in dss.keys():

                        dss[qcname] = ('time', 2 * np.ones(len(dss[name]), np.int8))
                        dss[qcname].attrs = utils.fill_required_qcattrs({}, name)
                        # 2 is "not eval"
                # outname = outdir + '/' + utils.get_file_id(dss) + '.nc'
                _log.info('Writing %s', outname)
                timeunits = 'seconds since 1970-01-01T00:00:00Z'
                timecalendar = 'gregorian'
                try:
                    del dss.profile_time.attrs['_FillValue']
                    del dss.profile_time.attrs['units']
                except KeyError:
                    pass
                dss.to_netcdf(outname, encoding={'time': {'units': timeunits,
                                                          'calendar': timecalendar,
                                                          'dtype': 'float64'},
                                                          'profile_time':
                                                         {'units': timeunits,
                                                         '_FillValue': -99999.0,
                                                         'dtype': 'float64'},
                }

                                                         )

                # add traj_strlen using bare ntcdf to make IOOS happy
                with netCDF4.Dataset(outname, 'r+') as nc:
                    nc.renameDimension('string%d' % trajlen, 'traj_strlen')

def make_gridfiles(inname, outdir, deploymentyaml, *, fnamesuffix='', dz=1, starttime='1970-01-01'):
    """
    Turn a timeseries netCDF file into a vertically gridded netCDF.

    Parameters
    ----------
    inname : str or Path
        netcdf file to break into profiles

    outdir : str or Path
        directory to place profiles

    deploymentyaml : str or Path
        location of deployment yaml file for the netCDF file.  This should
        be the same yaml file that was used to make the timeseries file.

    dz : float, default = 1
        Vertical grid spacing in meters.

    Returns
    -------
    outname : str
        Name of gridded netCDF file. The gridded netCDF file has coordinates of
        'depth' and 'profile', so each variable is gridded in depth bins and by
        profile number.  Each profile has a time, latitude, and longitude.
    """
    try:
        os.mkdir(outdir)
    except FileExistsError:
        pass

    deployment = utils._get_deployment(deploymentyaml)

    profile_meta = deployment['profile_variables']

    ds = xr.open_dataset(inname, decode_times=True)
    ds = ds.where(ds.time > np.datetime64(starttime), drop=True)
    _log.info(f'Working on: {inname}')
    _log.debug(str(ds))
    _log.debug(str(ds.time[0]))
    _log.debug(str(ds.time[-1]))

    profiles = np.unique(ds.profile_index)
    profiles = [p for p in profiles if (~np.isnan(p) and not (p % 1)
                                        and (p > 0))]
    profile_bins = np.hstack((np.array(profiles) - 0.5, [profiles[-1]+0.5]))
    _log.debug(profile_bins)
    Nprofiles = len(profiles)
    _log.info(f'Nprofiles {Nprofiles}')
    depth_bins = np.arange(0, 1100.1, dz)
    depths = depth_bins[:-1] + 0.5
    xdimname = 'time'
    dsout = xr.Dataset(
        coords={'depth': ('depth', depths),
                'profile': (xdimname, profiles)})
    dsout['depth'].attrs = {'units': 'm',
                            'long_name': 'Depth',
                            'standard_name': 'depth',
                            'positive': 'down',
                            'coverage_content_type': 'coordinate',
                            'comment': 'center of depth bins'}

    ds['time_1970'] = ds.temperature.copy()
    ds['time_1970'].values = ds.time.values.astype(np.float64)
    for td in ('time_1970', 'longitude', 'latitude'):
        good = np.where(~np.isnan(ds[td]) & (ds['profile_index'] % 1 == 0))[0]
        dat, xedges, binnumber = stats.binned_statistic(
                ds['profile_index'].values[good],
                ds[td].values[good], statistic='mean',
                bins=[profile_bins])
        if td == 'time_1970':
            td = 'time'
            dat = dat.astype('timedelta64[ns]') + np.datetime64('1970-01-01T00:00:00')
        _log.info(f'{td} {len(dat)}')
        dsout[td] = (('time'), dat, ds[td].attrs)
    ds.drop('time_1970')
    good = np.where(~np.isnan(ds['time']) & (ds['profile_index'] % 1 == 0))[0]
    _log.info(f'Done times! {len(dat)}')
    dsout['profile_time_start'] = (
        (xdimname), dat, profile_meta['profile_time_start'])
    dsout['profile_time_end'] = (
        (xdimname), dat, profile_meta['profile_time_end'])

    for k in ds.keys():
        if k in ['time', 'profile', 'longitude', 'latitude', 'depth'] or 'time' in k:
            continue
        _log.info('Gridding %s', k)
        good = np.where(~np.isnan(ds[k]) & (ds['profile_index'] % 1 == 0))[0]
        if len(good) <= 0:
            continue
        if "average_method" in ds[k].attrs:
            average_method = ds[k].attrs["average_method"]
            ds[k].attrs["processing"] = (
                f"Using average method {average_method} for "
                f"variable {k} following deployment yaml.")
            if average_method == "geometric mean":
                average_method = stats.gmean
                ds[k].attrs["processing"] += (" Using geometric mean implementation "
                                              "scipy.stats.gmean")
        else:
            average_method = "mean"
        dat, xedges, yedges, binnumber = stats.binned_statistic_2d(
                ds['profile_index'].values[good],
                ds['depth'].values[good],
                values=ds[k].values[good], statistic=average_method,
                bins=[profile_bins, depth_bins])

        _log.debug(f'dat{np.shape(dat)}')
        dsout[k] = (('depth', xdimname), dat.T, ds[k].attrs)

        # fill gaps in data:
        dsout[k].values = utils.gappy_fill_vertical(dsout[k].values)

    # fix u and v, because they should really not be gridded...
    if (('water_velocity_eastward' in dsout.keys()) and
            ('u' in profile_meta.keys())):
        _log.debug(str(ds.water_velocity_eastward))
        dsout['u'] = dsout.water_velocity_eastward.mean(axis=0)
        dsout['u'].attrs = profile_meta['u']
        dsout['v'] = dsout.water_velocity_northward.mean(axis=0)
        dsout['v'].attrs = profile_meta['v']
        dsout = dsout.drop(['water_velocity_eastward',
                            'water_velocity_northward'])
    dsout.attrs = ds.attrs
    dsout.attrs.pop('cdm_data_type')
    # fix to be ISO parsable:
    if len(dsout.attrs['deployment_start']) > 18:
        dsout.attrs['deployment_start'] = dsout.attrs['deployment_start'][:19]
        dsout.attrs['deployment_end'] = dsout.attrs['deployment_end'][:19]
        dsout.attrs['time_coverage_start'] = dsout.attrs['time_coverage_start'][:19]
        dsout.attrs['time_coverage_end'] = dsout.attrs['time_coverage_end'][:19]
    # fix standard_name so they don't overlap!
    try:
        dsout['waypoint_latitude'].attrs.pop('standard_name')
        dsout['waypoint_longitude'].attrs.pop('standard_name')
        dsout['profile_time_start'].attrs.pop('standard_name')
        dsout['profile_time_end'].attrs.pop('standard_name')
    except:
        pass
    # set some attributes for cf guidance
    # see H.6.2. Profiles along a single trajectory
    # https://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/build/aphs06.html
    dsout.attrs['featureType'] = 'trajectoryProfile'
    dsout['profile'].attrs['cf_role'] = 'profile_id'
    dsout['mission_number'] = int(1)
    dsout['mission_number'].attrs['cf_role'] = 'trajectory_id'
    dsout = dsout.set_coords(['latitude', 'longitude', 'time'])
    for k in dsout:
        if k in ['profile', 'depth', 'latitude', 'longitude', 'time', 'mission_number']:
            dsout[k].attrs['coverage_content_type'] = 'coordinate'
        else:
            dsout[k].attrs['coverage_content_type'] = 'physicalMeasurement'


    outname = outdir + '/' + ds.attrs['deployment_name'] + '_grid' + fnamesuffix + '.nc'
    _log.info('Writing %s', outname)
    # timeunits = 'nanoseconds since 1970-01-01T00:00:00Z'
    dsout.to_netcdf(
        outname,
        encoding={'time': {'units': 'seconds since 1970-01-01T00:00:00Z',
                           '_FillValue': np.nan,
                           'calendar': 'gregorian',
                           'dtype': 'float64'}})
    _log.info('Done gridding')

    return outname


# aliases
extract_L0timeseries_profiles = extract_timeseries_profiles
make_L0_gridfiles = make_gridfiles


__all__ = ['extract_timeseries_profiles', 'make_gridfiles']
