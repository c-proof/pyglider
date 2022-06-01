import logging
import xarray as xr
import numpy as np
import pyglider.utils as utils
import os
import yaml
import netCDF4
import scipy.stats as stats

_log = logging.getLogger(__name__)


def extract_timeseries_profiles(inname, outdir, deploymentyaml):
    """
    """
    try:
        os.mkdir(outdir)
    except FileExistsError:
        pass

    with open(deploymentyaml) as fin:
        deployment = yaml.safe_load(fin)
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
            if not os.path.exists(outname):
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

                dss['profile_id'] = np.array(p*1.0)
                dss['profile_id'].attrs = profile_meta['profile_id']
                dss['profile_time'] = dss.time.mean()
                dss['profile_time'].attrs = profile_meta['profile_time']
                dss['profile_lon'] = dss.longitude.mean()
                dss['profile_lon'].attrs = profile_meta['profile_lon']
                dss['profile_lat'] = dss.latitude.mean()
                dss['profile_lat'].attrs = profile_meta['profile_lat']

                dss['lat'] = dss['latitude']
                dss['lon'] = dss['longitude']
                dss['platform'] = np.NaN
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

                dss['lat_uv'] = np.NaN
                dss['lat_uv'].attrs = profile_meta['lat_uv']
                dss['lon_uv'] = np.NaN
                dss['lon_uv'].attrs = profile_meta['lon_uv']
                dss['time_uv'] = np.NaN
                dss['time_uv'].attrs = profile_meta['time_uv']

                dss['instrument_ctd'] = np.NaN
                dss['instrument_ctd'].attrs = profile_meta['instrument_ctd']

                dss.attrs['date_modified'] = str(np.datetime64('now')) + 'Z'

                # ancillary variables::
                to_fill = ['temperature', 'pressure', 'conductivity',
                           'salinity', 'density', 'lon', 'lat', 'depth']
                for name in to_fill:
                    dss[name].attrs['ancillary_variables'] = name + '_qc'

                # outname = outdir + '/' + utils.get_file_id(dss) + '.nc'
                _log.info('Writing %s', outname)
                if 'units' in dss.profile_time.attrs:
                    dss.profile_time.attrs.pop('units')
                timeunits = 'seconds since 1970-01-01T00:00:00Z'
                dss.to_netcdf(outname, encoding={'time': {'units': timeunits}})

                # add traj_strlen using bare ntcdf to make IOOS happy
                with netCDF4.Dataset(outname, 'r+') as nc:
                    nc.renameDimension('string%d' % trajlen, 'traj_strlen')


def make_gridfiles(inname, outdir, deploymentyaml, dz=1):
    """
    """
    try:
        os.mkdir(outdir)
    except FileExistsError:
        pass

    with open(deploymentyaml) as fin:
        deployment = yaml.safe_load(fin)
    profile_meta = deployment['profile_variables']

    ds = xr.open_dataset(inname)
    _log.info(f'Working on: {inname}')
    _log.debug(str(ds))
    _log.debug(str(ds.time[0]))
    _log.debug(str(ds.time[-1]))

    profiles = np.unique(ds.profile_index)
    profiles = [p for p in profiles if (~np.isnan(p) and not (p % 1)
                                        and (p > 0))]
    profile_bins = np.hstack((np.array(profiles) - 0.5, [profiles[-1]+0.5]))

    Nprofiles = len(profiles)
    _log.info(f'Nprofiles {Nprofiles}')
    depth_bins = np.arange(0, 1100.1, dz)
    depths = depth_bins[:-1] + 0.5

    dsout = xr.Dataset(
        coords={'depth': ('depth', depths),
                'profile': ('time', profiles)})
    ds['time_1970'] = ds.temperature.copy()
    ds['time_1970'].values = ds.time.values.astype(np.float64)/1e9
    for td in ('time_1970', 'longitude', 'latitude'):
        good = np.where(~np.isnan(ds[td]) & (ds['profile_index'] % 1 == 0))[0]
        dat, xedges, binnumber = stats.binned_statistic(
                ds['profile_index'].values[good],
                ds[td].values[good], statistic='mean',
                bins=[profile_bins])
        if td == 'time_1970':
            td = 'time'
            dat = dat.astype('timedelta64[s]') + np.datetime64('1970-01-01T00:00:00')
        _log.info(f'{td} {len(dat)}')
        dsout[td] = (('time'), dat, ds[td].attrs)
    ds.drop('time_1970')
    good = np.where(~np.isnan(ds['time']) & (ds['profile_index'] % 1 == 0))[0]
    _log.info(f'Done times! {len(dat)}')
    dsout['profile_time_start'] = (
        ('time'), dat, profile_meta['profile_time_start'])
    dsout['profile_time_end'] = (
        ('time'), dat, profile_meta['profile_time_end'])

    for k in ds.keys():
        if k in ['time', 'longitude', 'latitude', 'depth'] or 'time' in k:
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
        dsout[k] = (('depth', 'time'), dat.T, ds[k].attrs)

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

    outname = outdir + '/' + ds.attrs['deployment_name'] + '_grid.nc'
    _log.info('Writing %s', outname)
    timeunits = 'seconds since 1970-01-01T00:00:00Z'
    dsout.to_netcdf(outname, encoding={'time': {'units': timeunits}})
    _log.info('Done gridding')

    return outname


# aliases
extract_L0timeseries_profiles = extract_timeseries_profiles
make_L0_gridfiles = make_gridfiles


__all__ = ['extract_timeseries_profiles', 'make_gridfiles']
