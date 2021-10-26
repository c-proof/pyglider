# -*- coding: utf-8 -*-
import datetime
import glob
import itertools
import logging
from math import floor, fmod
import numpy as np
import os
import re
import datetime
import xarray as xr
import yaml
import pyglider.utils as utils
import pandas as pd


_log = logging.getLogger(__name__)


def _outputname(f, outdir):
    fnout = os.path.basename(f)
    fns = fnout.split('.')
    fns = fns[:5]
    print(fns)
    fns[4] = '%04d' % int(fns[4])
    fns[1] = '%04d' % int(fns[1])
    fnout = ''
    for ff in fns:
        fnout += ff.lower() + '.'
    filenum = int(fns[4])
    return outdir + fnout + 'nc', filenum


def _needsupdating(ftype, fin, fout):
    if ftype == 'pld1':
        fout = fout[:-3]+'_arod.nc'
    if not os.path.isfile(fout):
        return True
    return (os.path.getmtime(fin) >= os.path.getmtime(fout))

def _sort(ds):
    return ds.sortby('time')

def raw_to_rawnc(indir, outdir, deploymentyaml, incremental=True, min_samples_in_file=5):
    """
    Convert seaexplorer text files to raw netcdf files.

    Parameters
    ----------
    indir : str
        Directory with the raw files are kept.  Recommend naming this
        direectory "raw"

    outdir : str
        Directory to write the matching ``*.nc`` files. Recommend ``rawnc``.

    deploymentyaml : str
        YAML text file with deployment information for this glider.

    incremental : bool, optional
        If *True* (default), only netcdf files that are older than the
        binary files are re-parsed.

    min_samples_in_file : int
        Minimum number of samples in a raw file to trigger writing a netcdf file.
        Defaults to 5

    Returns
    -------
    status : bool
        *True* success.

    Notes
    -----

    This process can be slow for many files.

    """

    print(outdir)
    for ftype in ['gli', 'pld1']:
        for rawsub in ['raw', 'sub']:
            d = indir + f'*.{ftype}.{rawsub}.*'
            files = glob.glob(d)
            fnum = np.zeros(len(files))
            # these files don't sort properly, but we can sort them here.
            for n, f in enumerate(files):
                p = os.path.basename(f).split('.')
                fnum[n] = p[4]
            inds = np.argsort(fnum)
            files = [files[ind] for ind in inds]
            print(files)

            if len(files) < 0:
                raise FileNotFoundError('No raw files found in %s' % indir)

            try:
                os.mkdir(outdir)
            except FileExistsError:
                pass

            badfiles = []
            for ind, f in enumerate(files):
                # output name:
                fnout, filenum = _outputname(f, outdir)
                _log.info(f'{f} to {fnout}')
                if not incremental or _needsupdating(ftype, f, fnout):
                    _log.info(f'Doing: {f} to {fnout}')
                    out = pd.read_csv(f, header=0, delimiter=';',
                                            parse_dates=True, index_col=0,
                                            dayfirst=True)
                    with out.to_xarray() as outx:

                        key = list(outx.coords.keys())[0]
                        outx = outx.rename({key:'time'})
                        # dumb time down to seconds since 1970-01-01
                        outx['time'] = outx['time'].astype(np.float64)/1e9
                        outx['time'].attrs['units'] = (
                            'seconds since 1970-01-01T00:00:00Z')
                        outx['fnum'] = ('time',
                            int(filenum) * np.ones(len(outx['time'])))
                        if ftype == 'gli':
                            outx.to_netcdf(fnout[:-3]+'.nc', 'w')
                        else:
                            # split into instruments: this is hardcoded for now
                            gpctd = {'GPCTD_CONDUCTIVITY', 'GPCTD_TEMPERATURE',
                                     'GPCTD_PRESSURE', 'GPCTD_DOF'}
                            flbbcd = {'FLBBCD_CHL_COUNT', 'FLBBCD_CHL_SCALED',
                                      'FLBBCD_BB_700_COUNT',
                                      'FLBBCD_BB_700_SCALED',
                                      'FLBBCD_CDOM_COUNT',
                                      'FLBBCD_CDOM_SCALED'}
                            arod = {'AROD_FT_TEMP', 'AROD_FT_DO'}
                            nav = {'NAV_RESOURCE', 'NAV_LONGITUDE',
                                   'NAV_LATITUDE', 'NAV_DEPTH'}
                            pldGPCTD = outx.where(np.isfinite(
                                outx.GPCTD_TEMPERATURE), drop=True)
                            pldGPCTD = pldGPCTD.drop_vars(flbbcd)
                            pldGPCTD = pldGPCTD.drop_vars(arod)
                            if pldGPCTD.indexes["time"].size > min_samples_in_file:
                                pldGPCTD.to_netcdf(fnout[:-3] + '_gpctd.nc', 'w',
                                                   unlimited_dims=['time'])
                            else:
                                _log.warning('Number of GPCTD data points too small. Skipping nc write')
                                                                
                            pldGPCTD = outx.where(np.isfinite(
                                outx.FLBBCD_CHL_COUNT), drop=True)
                            pldGPCTD = pldGPCTD.drop_vars(gpctd)
                            pldGPCTD = pldGPCTD.drop_vars(arod)
                            if pldGPCTD.indexes["time"].size > min_samples_in_file:
                                pldGPCTD.to_netcdf(fnout[:-3] + '_flbbcd.nc', 'w',
                                                   unlimited_dims=['time'])
                            else:
                                _log.warning('Number of FLBBCD data points too small. Skipping nc write')

                            pldGPCTD = outx.where(np.isfinite(
                                outx.AROD_FT_DO), drop=True)
                            pldGPCTD = pldGPCTD.drop_vars(gpctd)
                            pldGPCTD = pldGPCTD.drop_vars(flbbcd)
                            if pldGPCTD.indexes["time"].size > min_samples_in_file:
                                pldGPCTD.to_netcdf(fnout[:-3] + '_arod.nc', 'w',
                                                   unlimited_dims=['time'])
                            else:
                                _log.warning('Number of AROD data points too small. Skipping nc write')

            if len(badfiles) > 0:
                _log.warning('Some files could not be parsed:')
                for fn in badfiles:
                    _log.warning('%s', fn)
            _log.info('All done!')

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
    print('id', id)
    outgli = outdir + '/' + id + '-rawgli.nc'
    outpld = outdir + '/' + id + '-' + kind + 'pld.nc'

    _log.info('Opening *.gli.sub.*.nc multi-file dataset from %s', indir)
    files = sorted(glob.glob(indir+'/*.gli.sub.*.nc'))
    with xr.open_dataset(files[0], decode_times=False) as gli:
        for fil in files[1:]:
            try:
                with xr.open_dataset(fil, decode_times=False) as gli2:
                    gli = xr.concat([gli, gli2], dim='time')
            except:
                pass
        gli.to_netcdf(outgli)
    _log.info(f'Done writing {outgli}')

    if kind == 'boo':
        _log.info('Opening *.pld.sub.*.nc multi-file dataset')
        files = sorted(glob.glob(indir+'/*.pld1.'+kind+'.*.nc'))
        with xr.open_dataset(files[0], decode_times=False) as pld:
            for fil in files[1:]:
                try:
                    with xr.open_dataset(fil, decode_times=False) as pld2:
                        pld = xr.concat([pld, pld2], dim='time')
                except:
                    pass
            _log.info('Writing ' + outpld)
            pld.to_netcdf(outpld)
    else:
        from dask.diagnostics import ProgressBar

        _log.info('Working on gpctd')
        print(f'{indir}/*.{kind}*gpctd.nc')
        try:
            os.remove(outpld[:-5]+'_gpctd.nc')
        except:
            pass
        with xr.open_mfdataset(f'{indir}/*.{kind}.*gpctd.nc',  decode_times=False, parallel=False, lock=False, preprocess=_sort) as pld:
            print(pld)
            print('Writing')
            delayed_obj = pld.to_netcdf(outpld[:-5]+'_gpctd.nc', 'w', unlimited_dims=['time'], compute=False)
            with ProgressBar():
                results = delayed_obj.compute()
        _log.info('Working on flbbcd')
        try:
            os.remove(outpld[:-5]+'_flbbcd.nc')
        except:
            pass
        with xr.open_mfdataset(f'{indir}/*.{kind}.*flbbcd.nc', decode_times=False, lock=False, preprocess=_sort) as pld:
            pld.to_netcdf(outpld[:-5]+'_flbbcd.nc', 'w', unlimited_dims=['time'])
        _log.info('Working on AROD')
        try:
            os.remove(outpld[:-5]+'_arod.nc')
        except:
            pass
        with xr.open_mfdataset(f'{indir}/*.{kind}.*arod.nc', decode_times=False, lock=False, preprocess=_sort) as pld:
            pld = pld.coarsen(time=8, boundary='trim').mean()
            pld.to_netcdf(outpld[:-5]+'_arod.nc', 'w', unlimited_dims=['time'])

    _log.info('Done merge_rawnc')

    return


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

def raw_to_L0timeseries(indir, outdir, deploymentyaml, kind='raw',
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
    gli = xr.open_dataset(indir + '/' + id + '-rawgli.nc', decode_times=False)
    ctd = xr.open_dataset(indir + '/' + id + '-' + kind+ 'p_gpctd.nc',
                          decode_times=False)
    arod = xr.open_dataset(indir + '/' + id + '-' + kind+ 'p_arod.nc',
                      decode_times=False)
    flb = xr.open_dataset(indir + '/' + id + '-' + kind+ 'p_flbbcd.nc',
                      decode_times=False)

    # build a new data set based on info in `deployment.`
    # We will use ebd.m_present_time as the interpolant if the
    # variabel is in dbd.

    ds = xr.Dataset()
    attr = {}
    name = 'time'
    for atts in ncvar[name].keys():
        if atts != 'coordinates':
            attr[atts] = ncvar[name][atts]

    # the ctd will be our timebase.  It oversamples the nav data, but
    # mildly undersamples the optics and oxygen....
    indctd = np.where(~np.isnan(ctd.GPCTD_TEMPERATURE))[0]

    print('TIME', ctd['time'])
    ds[name] = (('time'), ctd[name].values[indctd], attr)
    print(ds['time'])
    thenames = list(ncvar.keys())
    print(thenames)
    thenames.remove('time')

    for name in thenames:
        _log.info('interpolating ' + name)
        if not('method' in ncvar[name].keys()):
            # variables that are in the data set or can be interpolated from it
            if 'conversion' in ncvar[name].keys():
                convert = getattr(utils, ncvar[name]['conversion'])
            else:
                convert = utils._passthrough
            sensorname = ncvar[name]['source']
            if sensorname in ctd.keys():
                _log.debug('sensorname %s', sensorname)
                val = convert(ctd[sensorname])
                val = _interp_pld_to_pld(ctd, ds, val, indctd)
                ncvar['method'] = 'linear fill'
            elif sensorname in arod.keys():
                _log.debug('sensorname %s', sensorname)
                val = convert(arod[sensorname])
                val = _interp_pld_to_pld(arod, ds, val, indctd)
                ncvar['method'] = 'linear fill'
            elif sensorname in flb.keys():
                _log.debug('sensorname %s', sensorname)
                val = convert(flb[sensorname])
                val = _interp_pld_to_pld(flb, ds, val, indctd)
                ncvar['method'] = 'linear fill'
            else:
                val = gli[sensorname]
                #val = utils._zero_screen(val)
        #        val[val==0] = np.NaN
                val = convert(val)
                print('Gli', gli)
                val = _interp_gli_to_pld(gli, ds, val, indctd)

            # make the attributes:
            ncvar[name].pop('coordinates', None)
            attrs = ncvar[name]
            attrs = utils.fill_required_attrs(attrs)
            ds[name] = (('time'), val.data, attrs)

    # fix lon and lat to be linearly interpolated between fixes
    good = np.where(np.abs(np.diff(ds.longitude)) +
                    np.abs(np.diff(ds.latitude)) > 0)[0] + 1
    ds['longitude'].values = np.interp(ds.time,
        ds.time[good], ds.longitude[good])
    ds['latitude'].values = np.interp(ds.time,
        ds.time[good], ds.latitude[good])

    # some derived variables:
    ds = utils.get_glider_depth(ds)
    ds = utils.get_distance_over_ground(ds)
    #    ds = utils.get_profiles(ds)
    ds = utils.get_profiles_new(ds,
            filt_time=profile_filt_time, profile_min_time=profile_min_time)

    ds = utils.get_derived_eos_raw(ds)

    ds = ds.assign_coords(longitude=ds.longitude)
    ds = ds.assign_coords(latitude=ds.latitude)
    ds = ds.assign_coords(depth=ds.depth)
    #ds = ds._get_distance_over_ground(ds)

    ds = utils.fill_metadata(ds, deployment['metadata'], device_data)

    # somehow this comes out unsorted:
    ds = ds.sortby(ds.time)

    start = ((ds['time'].values[0]).astype('timedelta64[s]') +
        np.datetime64('1970-01-01T00:00:00'))
    end = ((ds['time'].values[-1]).astype('timedelta64[s]')  +
        np.datetime64('1970-01-01T00:00:00'))

    ds.attrs['deployment_start'] = str(start)
    ds.attrs['deployment_end'] = str(end)

    try:
        os.mkdir(outdir)
    except:
        pass
    id0 = ds.attrs['deployment_name']
    outname = outdir + id0 +  '.nc'
    _log.info('writing %s', outname)
    ds.to_netcdf(outname, 'w')

    return outname

# alias:
raw_to_L1timeseries = raw_to_L0timeseries
