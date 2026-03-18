"""
Routines that are used for common processing of netcdf files after they have
been converted to standard timeseries.
"""

import logging
import os

import netCDF4
import numpy as np
import scipy.stats as stats
from scipy import signal 
import xarray as xr
import gsw  
import yaml
from datetime import date



import pyglider.utils as utils

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
        profiles = [p for p in profiles if (~np.isnan(p) and not (p % 1) and (p > 0))]
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
                    'deployment of a glider and may span multiple data files.'
                )
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
                dss['profile_id'].attrs['valid_min'] = np.int32(
                    dss['profile_id'].attrs['valid_min']
                )
                dss['profile_id'].attrs['valid_max'] = np.int32(
                    dss['profile_id'].attrs['valid_max']
                )

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
                comment = meta['glider_model'] + ' operated by ' + meta['institution']
                dss['platform'].attrs['comment'] = comment
                dss['platform'].attrs['id'] = (
                    meta['glider_name'] + meta['glider_serial']
                )
                dss['platform'].attrs['instrument'] = 'instrument_ctd'
                dss['platform'].attrs['long_name'] = (
                    meta['glider_model'] + dss['platform'].attrs['id']
                )
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
                to_fill = [
                    'temperature',
                    'pressure',
                    'conductivity',
                    'salinity',
                    'density',
                    'lon',
                    'lat',
                    'depth',
                ]
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
                dss.to_netcdf(
                    outname,
                    encoding={
                        'time': {
                            'units': timeunits,
                            'calendar': timecalendar,
                            'dtype': 'float64',
                        },
                        'profile_time': {
                            'units': timeunits,
                            '_FillValue': -99999.0,
                            'dtype': 'float64',
                        },
                    },
                )

                # add traj_strlen using bare ntcdf to make IOOS happy
                with netCDF4.Dataset(outname, 'r+') as nc:
                    nc.renameDimension('string%d' % trajlen, 'traj_strlen')

def make_gridfiles(inname,outdir,deploymentyaml,*,fnamesuffix='',depth_bins=None,dz=1,starttime='1970-01-01',maskfunction=None,interp_variables=None):
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

    maskfunction : callable or None, optional
        Function applied to the dataset before gridding, usually to choose what data will be set to NaN based on quality flags. 

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
    
    if maskfunction is not None:
        ds = maskfunction(ds)

    ds = ds.where(ds.time > np.datetime64(starttime), drop=True)
    _log.info(f'Working on: {inname}')
    _log.debug(str(ds))
    _log.debug(str(ds.time[0]))
    _log.debug(str(ds.time[-1]))

    profiles = np.unique(ds.profile_index)
    profiles = [p for p in profiles if (~np.isnan(p) and not (p % 1) and (p > 0))]
    profile_bins = np.hstack((np.array(profiles) - 0.5, [profiles[-1] + 0.5]))
    _log.debug(profile_bins)
    Nprofiles = len(profiles)
    _log.info(f'Nprofiles {Nprofiles}')

    if depth_bins is None:
        # calculate depth bins using dz
        depth_bins = np.arange(0, 1100.1, dz)
    else:
        # sanity check user-provided bins
        if (
            depth_bins.ndim != 1 
            or not np.all(np.isfinite(depth_bins))
            or not np.issubdtype(depth_bins.dtype, np.number)
        ):
            raise ValueError('Depth bins must be a 1D array of finite numbers')
        if len(depth_bins) < 2:
            raise ValueError('There must be at least two depth bins edges')
        if not np.all(np.diff(depth_bins) > 0):
            raise ValueError('Depth bin edges must be strictly increasing and non-overlapping')
    
    # calculate bin centers
    depths = 0.5*(depth_bins[:-1] + depth_bins[1:])
    _log.debug(f'depth bins and centers {depth_bins} {{depths}}')

    xdimname = 'time'
    dsout = xr.Dataset(
        coords={'depth': ('depth', depths), 'profile': (xdimname, profiles)}
    )
    dsout['depth'].attrs = {
        'units': 'm',
        'long_name': 'Depth',
        'standard_name': 'depth',
        'positive': 'down',
        'source': ds.depth.attrs["source"], 
        'coverage_content_type': 'coordinate',
        'comment': 'center of depth bins',
    }

    # Bin by profile index, for the mean time, lat, and lon values for each profile
    ds['time_1970'] = ds.temperature.copy()
    ds['time_1970'].values = ds.time.values.astype(np.float64)
    for td in ('time_1970', 'longitude', 'latitude'):
        good = np.where(~np.isnan(ds[td]) & (ds['profile_index'] % 1 == 0))[0]
        dat, xedges, binnumber = stats.binned_statistic(
            ds['profile_index'].values[good],
            ds[td].values[good],
            statistic='mean',
            bins=[profile_bins],
        )
        if td == 'time_1970':
            td = 'time'
            dat = dat.astype('timedelta64[ns]') + np.datetime64('1970-01-01T00:00:00')
        _log.info(f'{td} {len(dat)}')
        dsout[td] = (('time'), dat, ds[td].attrs)

    # Bin by profile index, for the profile start (min) and end (max) times
    profile_lookup = {'profile_time_start': "min", 'profile_time_end': "max"}
    good = np.where(~np.isnan(ds['time']) & (ds['profile_index'] % 1 == 0))[0]
    for td, bin_stat in profile_lookup.items():
        _log.debug(f'td, bin_stat {td}, {bin_stat}')
        dat, xedges, binnumber = stats.binned_statistic(
            ds['profile_index'].values[good],
            ds['time_1970'].values[good],
            statistic=bin_stat,
            bins=[profile_bins],
        )
        dat = dat.astype('timedelta64[ns]') + np.datetime64('1970-01-01T00:00:00')
        _log.info(f'{td} {len(dat)}')
        dsout[td] = ((xdimname), dat, profile_meta[td])

    ds = drop_vars('time_1970')
    _log.info(f'Done times!')

    for k in ds.keys():
        if k in ['time', 'profile', 'longitude', 'latitude', 'depth'] or 'time' in k:
            continue
        _log.info('Gridding %s', k)
        good = np.where(~np.isnan(ds[k]) & (ds['profile_index'] % 1 == 0))[0]
        
        if len(good) <= 0:
            continue        
        if 'QC_protocol' in ds[k].attrs.values():
            # QC variables are treated as discrete flags rather than continuous data.
            # If a variable has a QC_protocol attribute, it is gridded using the
            # maximum flag in each bin (e.g. any QC3 in a bin makes the gridded bin QC3).
            method = np.nanmax
        else:
            # variables are treated as continuous data.
            # If a variable has a average_method attribute, it is gridded using the
            # mean in each bin  
            if 'average_method' in ds[k].attrs.values():
                method = ds[k].attrs['average_method']
                if method == 'geometric mean':
                    method = stats.gmean
            else:
                method = 'mean'
        
        dat, xedges, yedges, binnumber = stats.binned_statistic_2d(
                ds['profile_index'].values[good],
                ds['depth'].values[good],
                values=ds[k].values[good],
                statistic=method,
                bins=[profile_bins, depth_bins],
            )

        _log.debug(f'dat{np.shape(dat)}')
        dsout[k] = (('depth', xdimname), dat.T, ds[k].attrs)

        if interp_variables is not None:
            dsout[k] = interp_variables(dsout[k],ds[k])

    # fix u and v, because they should really not be gridded...
    if ('water_velocity_eastward' in dsout.keys()) and ('u' in profile_meta.keys()):
        _log.debug(str(ds.water_velocity_eastward))
        dsout['u'] = dsout.water_velocity_eastward.mean(axis=0)
        dsout['u'].attrs = profile_meta['u']
        dsout['v'] = dsout.water_velocity_northward.mean(axis=0)
        dsout['v'].attrs = profile_meta['v']
        dsout = dsout.drop(['water_velocity_eastward', 'water_velocity_northward'])
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
    # remove, so they can be encoded later:
    try:
        dsout['profile_time_start'].attrs.pop('units')
        dsout['profile_time_end'].attrs.pop('units')
        dsout['profile_time_start'].attrs.pop('_FillValue')
        dsout['profile_time_end'].attrs.pop('_FillValue')
    except:
        pass

    # set some attributes for cf guidance
    # see H.6.2. Profiles along a single trajectory
    # https://cfconventions.org/Data/cf-conventions/cf-conventions-1.7/build/aphs06.html
    dsout.attrs['featureType'] = 'trajectoryProfile'
    dsout['profile'].attrs['cf_role'] = 'profile_id'
    dsout['mission_number'] = np.int32(1)
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
    time_encoding = {
        'units': 'seconds since 1970-01-01T00:00:00Z',
        '_FillValue': np.nan,
        'calendar': 'gregorian',
        'dtype': 'float64',
    }
    dsout.to_netcdf(
        outname,
        encoding={
            'time': time_encoding, 
            'profile_time_start': time_encoding, 
            'profile_time_end': time_encoding, 
        },
    )
    _log.info('Done gridding')

    return outname

def CPROOF_mask(ds):
    """Mask QC4 samples in data variables (set to NaN) so gridding ignores them.
    Does NOT shorten arrays or overwrite QC variables.
    """
    _log = logging.getLogger(__name__)

    ds = ds.copy() 

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
    
    # QC variables: fill interpolatable NaN gaps with 1
    if 'QC_protocol' in attr.attrs.values():
        interp = var.interpolate_na(dim="depth", method="nearest", max_gap=50)
        filled = np.isnan(var) & np.isfinite(interp)
        return xr.where(filled, 1, var)

    # Continuous variables: linear interpolation
    return var.interpolate_na(dim="depth", method="linear", max_gap=50)

def read_ctd_constants_from_yaml(deploymentyaml):
    """
    Read CTD correction constants from deployment YAML.
    Missing values are returned as None.
    """
    with open(deploymentyaml, "r") as f:
        data = yaml.safe_load(f) or {}
    ctd = data.get("glider_devices", {}).get("ctd", {})
    thermal = ctd.get("Thermal_lag_constants_[alpha,tau]", None)
    alpha = None
    tau = None
    if thermal is not None:
        if len(thermal) > 0:
            alpha = thermal[0]
        if len(thermal) > 1:
            tau = thermal[1]
    dTdC = ctd.get("dTdC", None)
    return {
        "alpha": alpha,
        "tau": tau,
        "dTdC": dTdC,
    }
def ctd_constants(deployfile, *, alpha=None, tau=None, dTdC=None):
    """
    Read CTD correction constants from YAML, allow kwargs to override,
    warn on conflicts, and error if any required value is still missing.
    """
    with open(deployfile, "r") as file:
        data = yaml.safe_load(file) or {}

    atr = data.get("glider_devices", {}).get("ctd", {})

    thermal = atr.get("Thermal_lag_constants_[alpha,tau]")
    yaml_alpha = None
    yaml_tau = None
    if thermal is not None:
        if len(thermal) > 0:
            yaml_alpha = thermal[0]
        if len(thermal) > 1:
            yaml_tau = thermal[1]

    yaml_dTdC = atr.get("dTdC")

    vals_yaml = {"alpha": yaml_alpha, "tau": yaml_tau, "dTdC": yaml_dTdC}
    vals_kw = {"alpha": alpha, "tau": tau, "dTdC": dTdC}

    out = {}
    for key in vals_yaml:
        yaml_val = vals_yaml[key]
        kw_val = vals_kw[key]

        if kw_val is not None:
            if yaml_val is not None and yaml_val != kw_val:
                _log.warning(
                    "%s differs between YAML (%r) and kwargs (%r); using kwargs.",
                    key, yaml_val, kw_val
                )
            out[key] = kw_val
        else:
            out[key] = yaml_val

    missing = [k for k, v in out.items() if v is None]
    if missing:
        raise ValueError(
            f"Missing required CTD constants after checking kwargs and YAML: {missing}"
        )

    return out["alpha"], out["tau"], out["dTdC"]

def get_conductivity_clean(ts0, dT, dz, flag_stdev, clean_stdev, accuracy=None):

    """
    Temporarily flags any data points that are more than 5 standard deviations away from the overall time series mean 
    for a given depth bin and profile bin, then recomputes the mean and standard deviation, excluding the temporarily flagged values.
    Conductivity values that still differ from the mean by more than 3 standard deviations are flagged as 'bad' (QC 4). 
    If the difference between the 'bad' values and the mean is less than the accuracy of the sensor, then those points are not excluded.
    
    Parameters
    ----------
    ts0 : timeseries DataArray 
        Timeseries of mission  

    dT : float
        Ssize of the profile bins

    dz : float
        Size of depth bins in meters

    flag_stdev : float
        Number of standard deviations to temporarily flag bad salinity values 
        
    clean_stdev : float
        Number of standard deviations to flag bad conductivity values, after removing the temporary bad values from the calculation

    accuracy : callable or None, optional
        Accuracy of the sensor 
    
    Returns
    -------
    ts : timeseries 
        Timeseries of mission with new variable '
    """
    ts = ts0.copy(deep=True).load()
    ts = ts.where(np.isfinite(ts.conductivity), drop=False)
    Tbins = np.arange(np.min(ts.profile_index), np.max(ts.profile_index) + dT, dT)
    zbins = np.arange(np.min(ts.depth), np.max(ts.depth) + dz, dz)
    qc = np.full(ts.conductivity.shape, 4, dtype=int)
    finite_cond = np.isfinite(ts.conductivity.values)
    for n in range(len(Tbins) - 1):
        ind_Tbin = ((ts.profile_index.values >= Tbins[n]) &
            (ts.profile_index.values <= Tbins[n + 1]) &finite_cond)
        if not np.any(ind_Tbin):
            continue
        cond = ts.conductivity.values[ind_Tbin]
        depth = ts.depth.values[ind_Tbin]
        ind_bad_z = np.zeros(len(cond), dtype=bool)
        for m in range(len(zbins) - 1):
            ind_zbin = (depth >= zbins[m]) & (depth <= zbins[m + 1])
            if not np.any(ind_zbin):
                continue
            cond_z = cond[ind_zbin]
            cond_mean = np.nanmean(cond_z)
            cond_std = np.nanstd(cond_z)
            ind_flag = ((np.abs(cond_z - cond_mean) > flag_stdev * cond_std) &
                (np.abs(cond_z - cond_mean) > accuracy))
            # If everything got flagged, skip second-pass cleaning
            if np.all(ind_flag):
                ind_bad = ind_flag
            else:
                clean_mean = np.nanmean(cond_z[~ind_flag])
                clean_std = np.nanstd(cond_z[~ind_flag])
                ind_bad = ( (np.abs(cond_z - clean_mean) > clean_stdev * clean_std) &
                    (np.abs(cond_z - clean_mean) > accuracy))
            ind_bad_z[ind_zbin] = ind_bad
        qc_subset = np.where(ind_bad_z, 4, 1)# Good = 1, Bad = 4
        qc[ind_Tbin] = qc_subset
    ts = ts.assign(conductivity_QC=(ts.conductivity.dims, qc))
            
    return ts 


def CPROOF_sal_interpolate_filter(ds): 
    """
    Function applied to the dataset before finding the internal temperature. Function interpolates temperature over bad data and small data gaps
    to prevent errors from affecting the neighbouring cells. 

    Parameters
    ----------
    ds: DataArray
        Timeseries of mission data 
    
    Returns
    ----------
    interp: DataArray
        Timeseries of interpolated temperature
    
    """
    interp = ds["temperature"].where(ds["temperature_QC"] != 4)
    qc4 = (ds["temperature_QC"] == 4)
    qc4_buf = qc4.rolling(time=5, center=True, min_periods=1).max().astype(bool)
    interp = interp.where(~qc4_buf)

    interp = interp.interpolate_na(
        dim="time",
        method="linear",
        max_gap=np.timedelta64(60, "s"))

    return interp 
    
def correct_sal(ds, fn, alpha, tau, interpolate_filter = None):
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
        to prevent errors from affecting the neighbouring cells. 
    
    Returns 
    ----------
    sal: DataArray
        Timeseries of salinity_adjusted calculated using the internal temperature of the conductivity cell.
    """
    if interpolate_filter is not None: 
        temp = interpolate_filter(ds)
    else: 
        temp = ds.temperature_adjusted

    a = 4 * fn * alpha * tau / (1 + 4*fn*tau)
    b = 1 - 2 * a / alpha
    aa = [1, b]
    bb = [a, -a]
    tempcorr = temp.values.copy()
    tempcell = temp.values.copy()
    good = ~np.isnan(tempcell)
    tempcorr[good] = signal.lfilter(bb, aa, temp.values[good])
    tempcell = tempcell - tempcorr
    sal = gsw.SP_from_C(ds.conductivity* 10, tempcell, ds.pressure)

    return sal

def adjust_CTD(uncorrected_file,deploymentyaml,ts_directory, ds_directory, alpha= None, tau = None, dTdC=None, accuracy = None): 

    ts = xr.open_dataset(uncorrected_file)
    
    #############Identify anomalous conductivity values ############
    flag_stdev = 5 #number of standard deviations to temporarily flag bad salinity values 
    clean_stdev = 3 #number of standard deviations to flag bad conductivity values, after removing the temporary bad values from the calc
    dT = 50 #size of the profile bins
    dz = 5 #size of the depth bins

    ts0 = ts.copy() 
    ts0.conductivity[ts0.conductivity<0.1] = np.nan
    ts = get_conductivity_clean(ts0, dT, dz, flag_stdev, clean_stdev, accuracy=0.0003)

    ############## Flag salinity as QC4 where there are QC4 conductivity values (bad in == bad out) 
    sal_QC = xr.where(ts.conductivity_QC == 4,4,1) 
    ts['salinity_QC'] = sal_QC
    ts['temperature_QC'] =(('time'), np.ones_like(ts.temperature) )

    #added to C-PROOF files for gridding function 
    vars_ = ['conductivity_QC','salinity_QC','temperature_QC']
    for var in vars_:
        ts[var].attrs["average_method"] = "QC_protocol"

    ############## Resolve constants from kwargs and YAML
    alpha, tau, dTdC = ctd_constants(deploymentyaml,alpha=alpha,tau=tau,dTdC=dTdC)

    temp_adj = ts.temperature.copy()
    temp_adj.attrs = ts.temperature.attrs.copy()
    temp_adj.attrs['comment'] = 'temperature [degC]'

    ############## if there is a time lag 
    if dTdC is not None:
        dt = np.timedelta64(dTdC, 's')
        temp_adj = temp_adj.interp(time=ts.time + dt)
        temp_adj.attrs = ts.temperature.attrs.copy()
        temp_adj.attrs['comment'] = 'temperature [degC]'
        temp_adj.attrs['time_lag'] = f'{dTdC} second CT lag corrected'
        ts.attrs['dTdC'] = f'{dTdC} second CT lag corrected'

    ts['temperature_adjusted'] = temp_adj

    ############## Build salinity_adjusted
    if tau is not None:
        dt = np.diff(ts.time.values).astype("timedelta64[s]").astype(int)
        vals, counts = np.unique(dt, return_counts=True)
        srate = vals[np.argmax(counts)]
        fs = 1 / float(srate)
        fn = 0.5 * fs

        s = correct_sal(ts,fn,alpha=alpha,tau=tau,)

        sal_adj = xr.where(ts.salinity_QC == 1, s, ts.salinity)
        sal_adj.attrs = ts.salinity.attrs.copy()
        sal_adj.attrs['thermal_lag'] = (f'adjusted salinity [psu] using a thermal lag correction with alpha = {alpha} and tau = {tau}')
        if dTdC is not None:
            sal_adj.attrs['time_lag'] = f'found using temperature_adjusted corrected for {dTdC} second CT lag'

        ts['salinity_adjusted'] = sal_adj
        ts.attrs['correction_constants'] = f'alpha = {alpha}; tau={tau}'

    else:
        sal_adj = xr.DataArray(
            gsw.conversions.SP_from_C(10 *ts['conductivity'],ts['temperature_adjusted'],ts.pressure).values,dims=ts.salinity.dims,coords=ts.salinity.coords)
        sal_adj.attrs = ts.salinity.attrs.copy()
        if dTdC is not None:
            sal_adj.attrs['time_lag'] = f'found using temperature_adjusted'

        ts['salinity_adjusted'] = sal_adj
    
    ############# Update metadata 
    ts.attrs['processing_details'] = 'Processing details are located on the C-PROOF website for this mission under the reports tab.'  
    ts.attrs['processing_tech'] = 'Lauryn Talbot; ltalbot@uvic.ca'  
    ts.attrs['citation'] = '"Klymak, J., & Ross, T. (2025). C-PROOF Underwater Glider Deployment Datasets [Data set]. Canadian-Pacific Robotic Ocean Observing Facility.doi:10.82534/44DS-K310"'
    ts.attrs['references'] = 'https://doi.org/10.82534/44DS-K310' 
    ts.attrs['quality_flags']=['1 = good data; 3 = bad data, potentially correctable; 4 = bad data; 8 = estimated data']
    
    # Uncorrected conductivity
    ts['conductivity'].attrs['comment'] = 'uncorrected conductivity'
    ts['conductivity_QC'] = ts.conductivity_QC
        
    # Uncorrected temperature
    ts['temperature'].attrs['comment'] = 'uncorrected temperature [degC]'
    ts['temperature_QC'] = ts.temperature_QC 
    
    # Adjusted temperature
    ts['temperature_adjusted_QC'] =  ts['temperature_QC']  
    
    # Uncorrected salinity
    ts['salinity'].attrs['comment'] = 'uncorrected salinity [psu]'
    
    # Corrected salinity
    ts['salinity_adjusted_QC'] =  ts['salinity_QC'] 
        
    # Unadjusted density
    ts['density'].attrs['comment'] = 'unadjusted density'
    ts['density_QC'] =  ts['salinity_QC'] #same as salinity since derived 

    ########## Recalculate derived variables

    long = ts.longitude.fillna(ts.longitude.mean(skipna=True))
    lat = ts.latitude.fillna(ts.latitude.mean(skipna=True))
    sa_adj = gsw.SA_from_SP(ts['salinity_adjusted'],ts['pressure'],long,lat)
    ct_adj = gsw.CT_from_t(sa_adj,ts['temperature_adjusted'],ts['pressure'])
    ts['potential_density_adjusted'] = (('time'),1000 + gsw.density.sigma0(sa_adj,ct_adj).values)
    ts['potential_density_adjusted'].attrs['comment'] = 'calculated using adjusted salinity'
    ts['potential_density_adjusted_QC'] = ts['salinity_adjusted_QC']
    
    ts['potential_temperature_adjusted'] = (('time'),
                                            gsw.conversions.pt0_from_t(ts.salinity_adjusted,ts.temperature_adjusted,ts.pressure).values)
    ts['potential_temperature_adjusted'].attrs['comment'] = 'calculated using adjusted salinity'
    ts['potential_temperature_adjusted_QC'] = ts['salinity_adjusted_QC']

    processing_date = date.today().strftime('%Y%m%d')
    vars_ = ['salinity_adjusted','temperature_adjusted','potential_density_adjusted','potential_temperature_adjusted']
    for var in vars_:
        ts[var].attrs['processing_date'] = processing_date

    vars_ = ['conductivity_QC','salinity_QC','temperature_QC','salinity_adjusted_QC',
             'temperature_adjusted_QC','potential_density_adjusted_QC','potential_temperature_adjusted']
    for var in vars_:
        ts[var].attrs["comment"] = ['1 = good data; 3 = bad data, potentially correctable; 4 = bad data; 8 = estimated data']

    ########### Make gridfile 
    deploy_name = ts.deployment_name
    ts.to_netcdf(f'{ts_directory}/{deploy_name}_CTDadjusted.nc')
    make_gridfiles(f'{ts_directory}/{deploy_name}_CTDadjusted.nc',  f'{ds_directory}', deploymentyaml, fnamesuffix='_CTDadjusted',
                       maskfunction=CPROOF_mask,interp_variables=interpolate_vertical) 

    return 'Files are saved' 



# aliases
extract_L0timeseries_profiles = extract_timeseries_profiles
make_L0_gridfiles = make_gridfiles


__all__ = ['extract_timeseries_profiles', 'make_gridfiles']
