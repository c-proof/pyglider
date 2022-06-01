# -*- coding: utf-8 -*-
"""
SeaExplorer-specific processing routines.
"""
import glob
import logging
import numpy as np
import os
import xarray as xr
import yaml
import pyglider.utils as utils
import pandas as pd


_log = logging.getLogger(__name__)


def _outputname(f, outdir):
    fnout = os.path.basename(f)
    fns = fnout.split('.')
    fns = fns[:5]
    fns[4] = '%04d' % int(fns[4])
    fns[1] = '%04d' % int(fns[1])
    fnout = ''
    for ff in fns:
        fnout += ff.lower() + '.'
    filenum = int(fns[4])
    return outdir + fnout + 'nc', filenum


def _needsupdating(ftype, fin, fout):
    if not os.path.isfile(fout):
        return True
    return (os.path.getmtime(fin) >= os.path.getmtime(fout))


def _sort(ds):
    return ds.sortby('time')


def raw_to_rawnc(indir, outdir, deploymentyaml, incremental=True,
                 min_samples_in_file=5):
    """
    Convert seaexplorer text files to raw netcdf files.

    Parameters
    ----------
    indir : str
        Directory with the raw files are kept.  Recommend naming this
        directory "raw"

    outdir : str
        Directory to write the matching ``*.nc`` files. Recommend ``rawnc``.

    deploymentyaml : str
        YAML text file with deployment information for this glider.

    incremental : bool, optional
        If *True* (default), only netcdf files that are older than the
        binary files are re-parsed.

    min_samples_in_file : int
        Minimum number of samples in a raw file to trigger writing a netcdf
        file. Defaults to 5

    Returns
    -------
    status : bool
        *True* success.

    Notes
    -----

    This process can be slow for many files.

    """
    # Create out directory for netcdfs if it does not exist
    try:
        os.mkdir(outdir)
    except FileExistsError:
        pass

    for ftype in ['gli', 'pld1']:
        for rawsub in ['raw', 'sub']:
            _log.info(f'Reading in raw files matching *{ftype}.{rawsub}*')
            d = indir + f'*.{ftype}.{rawsub}.*'
            files = glob.glob(d)
            fnum = np.zeros(len(files))
            # these files don't sort properly, but we can sort them here.
            for n, f in enumerate(files):
                p = os.path.basename(f).split('.')
                fnum[n] = p[4]
            inds = np.argsort(fnum)
            files = [files[ind] for ind in inds]
            _log.info(f'found {files}')

            if len(files) == 0:
                # If no files of this type found, try the next type
                continue

            badfiles = []
            goodfiles = []
            for ind, f in enumerate(files):
                # output name:
                fnout, filenum = _outputname(f, outdir)
                _log.info(f'{f} to {fnout}')
                if not incremental or _needsupdating(ftype, f, fnout):
                    _log.info(f'Doing: {f} to {fnout}')
                    out = pd.read_csv(f, header=0, delimiter=';',
                                      parse_dates=True, index_col=0,
                                      dayfirst=True)
                    # If AD2CP data present, convert timestamps to datetime
                    if 'AD2CP_TIME' in out.columns:
                        # Set datestamps with date 00000 to None
                        out.loc[out.AD2CP_TIME.str[:6] == '000000',
                                'AD2CP_TIME'] = None
                        out['AD2CP_TIME'] = pd.to_datetime(out.AD2CP_TIME)
                    with out.to_xarray() as outx:
                        key = list(outx.coords.keys())[0]
                        outx = outx.rename({key: 'time'})
                        outx['fnum'] = ('time',
                                        int(filenum)*np.ones(len(outx['time'])))
                        if ftype == 'gli':
                            outx.to_netcdf(fnout[:-3] + '.nc', 'w')
                            goodfiles.append(f)
                        else:
                            if outx.indexes["time"].size > min_samples_in_file:
                                outx.to_netcdf(f'{fnout[:-3]}.nc', 'w',
                                               unlimited_dims=['time'])
                                goodfiles.append(f)
                            else:
                                _log.warning('Number of sensor data points'
                                             'too small. Skipping nc write')
                                badfiles.append(f)
            if len(badfiles) > 0:
                _log.warning('Some files could not be parsed:')
                for fn in badfiles:
                    _log.warning('%s', fn)
            if not goodfiles:
                _log.warning(f'No valid unprocessed seaexplorer files found in'
                             f'{indir}')
                return False
    _log.info('All raw files converted to nc')
    return True


def drop_rogue_1970(ds):
    """
    If dates greater than 1971, 1, 1 are observed, drop any dates before
    1971-01-01, from the datset and return it. This function removes 1970
    timestamps caused by a SeaExplorer rebooting during a mission. If all dates
    are < 1971-01-01, no action is taken

    Parameters:
        ds: xarray.DataSet
            dataset to check for pre-1971 dates
    Returns:
        ds: xarray.DataSet
    """
    dt_1971 = np.datetime64("1971-01-01")
    # If all dates before or after 1971-01-01, return the dataset
    if (ds.time > dt_1971).all() or (ds.time < dt_1971).all():
        return ds
    return ds.where(ds.time > dt_1971, drop=True)


def merge_rawnc(indir, outdir, deploymentyaml, incremental=False, kind='raw'):
    """
    Merge all the raw netcdf files in indir.  These are meant to be
    the raw flight and science files from the slocum.

    Parameters
    ----------
    indir : str
        Directory where the raw ``*.ebd.nc`` and ``*.dbd.nc`` files are.
        Recommend: ``./rawnc``

    outdir : str
        Directory where merged raw netcdf files will be put. Recommend:
        ``./rawnc/``.  Note that the netcdf files will be named following
        the data in *deploymentyaml*:
        ``glider_nameglider_serial-YYYYmmddTHHMM-rawebd.nc`` and
        ``...rawdbd.nc``.

    deploymentyaml : str
        YAML text file with deployment information for this glider.

    incremental : bool
        Only add new files....
    """

    with open(deploymentyaml) as fin:
        deployment = yaml.safe_load(fin)
    metadata = deployment['metadata']
    id = metadata['glider_name']
    outgli = outdir + '/' + id + '-rawgli.nc'
    outpld = outdir + '/' + id + '-' + kind + 'pld.nc'

    _log.info('Opening *.gli.sub.*.nc multi-file dataset from %s', indir)
    files = sorted(glob.glob(indir+'/*.gli.sub.*.nc'))
    if not files:
        _log.warning(f'No *gli*.nc files found in {indir}')
        return False
    with xr.open_dataset(files[0]) as gli:
        for fil in files[1:]:
            try:
                with xr.open_dataset(fil) as gli2:
                    gli = xr.concat([gli, gli2], dim='time')
            except:
                pass
        gli = drop_rogue_1970(gli)
        gli.to_netcdf(outgli)
    _log.info(f'Done writing {outgli}')

    _log.info('Opening *.pld.sub.*.nc multi-file dataset')
    files = sorted(glob.glob(indir+'/*.pld1.'+kind+'.*.nc'))
    if not files:
        _log.warning(f'No *{kind}*.nc files found in {indir}')
        return False
    with xr.open_dataset(files[0]) as pld:
        for fil in files[1:]:
            try:
                with xr.open_dataset(fil) as pld2:
                    pld = xr.concat([pld, pld2], dim='time')
            except:
                pass
        pld = drop_rogue_1970(pld)
        pld.to_netcdf(outpld)
    _log.info(f'Done writing {outpld}')
    _log.info('Done merge_rawnc')
    return True


def _interp_gli_to_pld(gli, ds, val, indctd):
    gli_ind = ~np.isnan(val)
    valout = np.interp(ds['time'],
                       gli['time'][gli_ind],
                       val[gli_ind])
    return valout


def _interp_pld_to_pld(pld, ds, val, indctd):
    pld_ind = np.where(~np.isnan(val))[0]
    if len(pld_ind) != len(indctd):
        val = np.interp(ds['time'],
                        pld['time'][pld_ind],
                        val[pld_ind])
    else:
        val = val[indctd]
    return val


def raw_to_timeseries(indir, outdir, deploymentyaml, kind='raw',
                      profile_filt_time=100, profile_min_time=300):
    """
    A little different than above, for the 4-file version of the data set.
    """

    with open(deploymentyaml) as fin:
        deployment = yaml.safe_load(fin)
    metadata = deployment['metadata']
    ncvar = deployment['netcdf_variables']
    device_data = deployment['glider_devices']
    id = metadata['glider_name']
    _log.info(f'Opening combined nav file {indir}/{id}-rawgli.nc')
    gli = xr.open_dataset(f'{indir}/{id}-rawgli.nc')
    _log.info(f'Opening combined payload file {indir}/{id}-{kind}pld.nc')
    sensor = xr.open_dataset(f'{indir}/{id}-{kind}pld.nc')

    # build a new data set based on info in `deploymentyaml.`
    # We will use ctd as the interpolant
    ds = xr.Dataset()
    attr = {}
    name = 'time'
    for atts in ncvar[name].keys():
        if atts != 'coordinates':
            attr[atts] = ncvar[name][atts]

    # If present, use the timebase specified in ncva: timebase in the
    # deployment yaml. Otherwise, the ctd will be our timebase.
    # It oversamples the nav data, but mildly undersamples the optics and
    # oxygen....
    if 'timebase' in ncvar:
        indctd = np.where(~np.isnan(sensor[ncvar['timebase']['source']]))[0]
    elif 'GPCTD_TEMPERATURE' in list(sensor.variables):
        _log.warning('No timebase specified. Using GPCTD_TEMPERATURE as time'
                     'base')
        indctd = np.where(~np.isnan(sensor.GPCTD_TEMPERATURE))[0]
    elif 'LEGATO_TEMPERATURE' in list(sensor.variables):
        _log.warning('No timebase specified. Using LEGATO_TEMPERATURE as time'
                     'base')
        indctd = np.where(~np.isnan(sensor.LEGATO_TEMPERATURE))[0]
    else:
        _log.warning('No gpctd or legato data found. Using NAV_DEPTH as time'
                     'base')
        indctd = np.where(~np.isnan(sensor.NAV_DEPTH))[0]
    ds['time'] = (('time'), sensor['time'].values[indctd], attr)
    thenames = list(ncvar.keys())
    for i in ['time', 'timebase', 'keep_variables']:
        if i in thenames:
            thenames.remove(i)
    for name in thenames:
        _log.info('interpolating ' + name)
        if not('method' in ncvar[name].keys()):
            # variables that are in the data set or can be interpolated from it
            if 'conversion' in ncvar[name].keys():
                convert = getattr(utils, ncvar[name]['conversion'])
            else:
                convert = utils._passthrough
            sensorname = ncvar[name]['source']
            if sensorname in list(sensor.variables):
                _log.debug('sensorname %s', sensorname)
                val = convert(sensor[sensorname])
                if 'coarsen' in ncvar[name]:
                    # smooth oxygen data as originally perscribed
                    coarsen_time = ncvar[name]['coarsen']
                    sensor_sub = sensor.coarsen(time=coarsen_time,
                                                boundary='trim').mean()
                    val2 = sensor_sub[sensorname]
                    val = _interp_gli_to_pld(sensor_sub, sensor, val2, indctd)
                val = val[indctd]

                ncvar['method'] = 'linear fill'
            else:
                val = gli[sensorname]
                val = convert(val)
                # Values from the glider netcdf must be interpolated to match
                # the sensor netcdf
                val = _interp_gli_to_pld(gli, ds, val, indctd)

            # make the attributes:
            ncvar[name].pop('coordinates', None)
            attrs = ncvar[name]
            attrs = utils.fill_required_attrs(attrs)
            ds[name] = (('time'), val.data, attrs)

    # fix lon and lat to be linearly interpolated between fixes
    good = np.where(np.abs(np.diff(ds.longitude)) +
                    np.abs(np.diff(ds.latitude)) > 0)[0] + 1
    ds['longitude'].values = np.interp(ds.time, ds.time[good],
                                       ds.longitude[good])
    ds['latitude'].values = np.interp(ds.time, ds.time[good],
                                      ds.latitude[good])

    # keep only timestamps with data from one of a set of variables
    if 'keep_variables' in ncvar:
        keeps = np.empty(len(ds.longitude))
        keeps[:] = np.nan
        keeper_vars = ncvar['keep_variables']
        for keep_var in keeper_vars:
            keeps[~np.isnan(ds[keep_var].values)] = 1
        ds = ds.where(~np.isnan(keeps))
        ds = ds.dropna(dim='time', how='all')

    # some derived variables:
    ds = utils.get_glider_depth(ds)
    ds = utils.get_distance_over_ground(ds)
    #    ds = utils.get_profiles(ds)
    ds = utils.get_profiles_new(ds, filt_time=profile_filt_time,
                                profile_min_time=profile_min_time)
    ds = utils.get_derived_eos_raw(ds)

    # Correct oxygen if present:
    if 'oxygen_concentration' in ncvar.keys():
        if 'correct_oxygen' in ncvar['oxygen_concentration'].keys():
            ds = utils.oxygen_concentration_correction(ds, ncvar)
        else:
            _log.warning('correct_oxygen not found in oxygen yaml. No'
                         'correction applied')
    ds = ds.assign_coords(longitude=ds.longitude)
    ds = ds.assign_coords(latitude=ds.latitude)
    ds = ds.assign_coords(depth=ds.depth)
    # ds = ds._get_distance_over_ground(ds)

    ds = utils.fill_metadata(ds, deployment['metadata'], device_data)

    # somehow this comes out unsorted:
    ds = ds.sortby(ds.time)

    start = ds['time'].values[0]
    end = ds['time'].values[-1]

    ds.attrs['deployment_start'] = str(start)
    ds.attrs['deployment_end'] = str(end)

    try:
        os.mkdir(outdir)
    except:
        pass
    id0 = ds.attrs['deployment_name']
    outname = outdir + id0 + '.nc'
    _log.info('writing %s', outname)
    if 'units' in ds.time.attrs.keys():
        ds.time.attrs.pop('units')
    if 'calendar' in ds.time.attrs.keys():
        ds.time.attrs.pop('calendar')
    if 'ad2cp_time' in list(ds):
        if 'units' in ds.ad2cp_time.attrs.keys():
            ds.ad2cp_time.attrs.pop('units')
    ds.to_netcdf(outname, 'w',
                 encoding={'time': {'units':
                                    'seconds since 1970-01-01T00:00:00Z'}})
    return outname


# alias:
raw_to_L1timeseries = raw_to_L0timeseries = raw_to_timeseries


def _parse_sensor_config(filename):
    """
    Reads the sensor config file of a SeaExplorer and extracts the active
    sensors and their calibration data.

    Parameters
    ----------
    filename: path to seapayload.cfg file

    Returns
    -------

    Dictionary of devices and their metadata as dictionaries of key_value pairs
    """
    # filename =
    # "/home/callum/Documents/data-flow/data-in/SEA063_M22/2_PLD/configs/seapayload.cfg"
    file = open(filename, 'r').read().split('\n')
    devices = []
    device_id = 'dummy_value'
    device_dicts = {}
    dict_for_device = {}
    for line in file:
        # Strip trailing whitespace
        line = line.strip(" ")
        # Look for key:value pairs
        if '=' in line and ">" not in line:
            # Split only on first = or AD2CP will break the parser
            key, value = line.split("=", 1)
            # Look for non-empty device declarations
            if key == "device" and value != "":
                devices.append(value)
            else:
                dict_for_device[key] = value
        # Look for pattern [devicename]
        elif line[1:-1] in devices:
            # add previous device to the device dict
            device_dicts[device_id] = dict_for_device
            device_id = line[1:-1]
            dict_for_device = {}
    # Append the final device to the dict
    device_dicts[device_id] = dict_for_device

    active_device_dicts = {k: device_dicts[k] for k in
                           (device_dicts.keys() & {*devices})}
    return active_device_dicts


__all__ = ['raw_to_rawnc', 'merge_rawnc', 'raw_to_timeseries']
