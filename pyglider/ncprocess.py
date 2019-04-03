import collections
import datetime
import seawater
import xarray as xr
import numpy as np
import pyglider.utils as utils
import os
import yaml
import netCDF4

def extract_L1timeseries_profiles(inname, outdir, deploymentyaml):
    """
    """
    try:
        os.mkdir(outdir)
    except FileExistsError:
        pass

    with open(deploymentyaml) as fin:
        deployment = yaml.safe_load(fin)
    meta = deployment['metadata']

    with xr.open_dataset(inname, decode_times=False) as ds:
        print(ds)
        profiles = np.unique(ds.profile_index)
        profiles = [p for p in profiles if (~np.isnan(p) and not (p % 1)
                                            and (p > 0))]
        for p in profiles:

            ind = np.where(ds.profile_index==p)[0]
            dss = ds.isel(time=ind)
            # this is the id for the whole file, not just this profile..
            dss['trajectory'] =  utils.get_file_id(ds).encode()
            trajlen = len(utils.get_file_id(ds).encode())
            dss['trajectory'].attrs['cf_role'] = 'trajectory_id'
            dss['trajectory'].attrs['comment'] = ('A trajectory is a single'
                    'deployment of a glider and may span multiple data files.')
            dss['trajectory'].attrs['long_name'] = \
                    'Trajectory/Deployment Name'

            # profile-averaged variables....
            profile_meta = deployment['profile_variables']
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
            comment = (meta['glider_model'] + ' opperated by ' +
                       meta['institution'])
            dss['platform'].attrs['comment'] =  comment
            dss['platform'].attrs['id'] =  (meta['glider_name'] +
                                            meta['glider_serial'])
            dss['platform'].attrs['instrument'] =  'instrument_ctd'
            dss['platform'].attrs['long_name'] = (meta['glider_model'] +
                    dss['platform'].attrs['id'])
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

            print('Hi', dss[name])


            outname = outdir + '/' + utils.get_file_id(dss) + '.nc'
            dss.to_netcdf(outname)
            # add traj_strlen using bare ntcdf..
            with netCDF4.Dataset(outname, 'r+') as nc:
                nc.renameDimension('string%d' % trajlen, 'traj_strlen')


def make_L2_gridfiles(inname, outdir, deploymentyaml):
    """
    """
