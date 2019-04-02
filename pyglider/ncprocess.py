import collections
import seawater
import xarray as xr
import numpy as np
import pyglider.utils as utils
import os
import yaml

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
            dss['trajectory'].attrs['cf_role'] = 'trajectory_id'
            dss['trajectory'].attrs['comments'] = ('A trajectory is a single'
                    'deployment of a glider and may span multiple data files.')
            dss['trajectiory'].attrs['long_name'] = \
                    'Trajectory/Deployment Name'

            # get some variables....
            dss['u'] = dss.water_velocity_eastward.mean()
            dss['v'] = dss.water_velocity_northward.mean()
            dss['profile_id'] = p
            dss['profile_time'] = dss.time.mean()
            dss['profile_lon'] = dss.longitude.mean()
            dss['profile_lat'] = dss.latitude.mean()
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
            dss.attrs['wmo_id'] = meta['wmo_id']

            dss['lat_uv'] = np.NaN
            dss['lon_uv'] = np.NaN

            outname = outdir + '/' + utils.get_file_id(dss) + '.nc'
            dss.to_netcdf(outname)
