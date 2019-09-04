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


def _needsupdating(fin, fout):
    if not os.path.isfile(fout):
        return True
    return (os.path.getmtime(fin) >= os.path.getmtime(fout))


def raw_to_rawnc(indir, outdir, deploymentyaml, incremental=True):
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
                if not incremental or _needsupdating(f, fnout):
                    _log.info(f'Doing: {f} to {fnout}')
                    if 1:
                        out = pd.read_csv(f, header=0, delimiter=';',
                                            parse_dates=True, index_col=0,
                                            dayfirst=True)
                        out = out.to_xarray()
                        key = list(out.coords.keys())[0]
                        out = out.rename({key:'time'})
                        # dumb time down to seconds since 1970-01-01
                        out['time'] = out['time'].astype(np.float64)/1e9
                        out['time'].attrs['units'] = (
                            'seconds since 1970-01-01T00:00:00Z')
                        out['fnum'] = ('time',
                            int(filenum) * np.ones(len(out['time'])))
                        out.to_netcdf(fnout, 'w')
                    else:
                        badfiles += [files[ind]]
                        _log.warning('Could not do parsing for %s', files[ind])

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
    id = metadata['glider_name'] + metadata['glider_serial']
    print('id', id)
    outgli = outdir + '/' + id + '-rawgli.nc'
    outpld = outdir + '/' + id + '-' + kind + 'pld.nc'

    _log.info('Opening *.gli.sub.*.nc multi-file dataset from %s', indir)
    gli = xr.open_mfdataset(indir + '/*.gli.sub.*.nc', decode_times=False)

    _log.info('Opening *.pld.sub.*.nc multi-file dataset')
    pld = xr.open_mfdataset(indir + '/*.pld1.' + kind + '.*.nc', decode_times=False)

    _log.info('Opening existing merged *.ebd.nc')

    # don't write if output files are up to date:
    # maybe this would be faster w/ a merge?  Not sure.
    write = False
    try:
        with xr.open_dataset(outgli, decode_times=False) as biggli:
            ind = np.where(gli.time > biggli.time[-1])[0]
            if len(ind) > 0:
                _log.info('New data found in raw *.ebd.nc')
                write = True
    except FileNotFoundError:
        _log.info('Merged file not found: %s', outgli )
        write = True
    if write:
        _log.info('Writing ' + outgli)
        gli.to_netcdf(outgli)
    else:
        _log.info('Not overwriting ' + outgli)
    _log.info('Done writing ')
    try:
        with xr.open_dataset(outpld, decode_times=False) as bigpld:
            ind = np.where(pld.time > bigpld.time[-1])[0]
            if len(ind) > 0:
                _log.info('New data found in raw *.dbd.nc')
                write = True
    except FileNotFoundError:
        _log.info('Merged file not found: %s', outpld )
        write = True
    if write:
        _log.info('Writing ' + outpld)
        pld.to_netcdf(outpld)
    else:
        _log.info('Not overwriting ' + outpld)

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


def raw_to_L1timeseries(indir, outdir, deploymentyaml, kind='raw'):
    """
    """

    with open(deploymentyaml) as fin:
        deployment = yaml.safe_load(fin)
    metadata = deployment['metadata']
    ncvar = deployment['netcdf_variables']

    id = metadata['glider_name'] + metadata['glider_serial']
    gli = xr.open_dataset(indir + '/' + id + '-rawgli.nc', decode_times=False)
    pld = xr.open_dataset(indir + '/' + id + '-' + kind+ 'pld.nc',
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
    indctd = np.where(~np.isnan(pld.GPCTD_TEMPERATURE))[0]

    ds[name] = (('time'), pld[name].values[indctd], attr)

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
            if sensorname in pld.keys():
                _log.debug('sensorname %s', sensorname)
                val = convert(pld[sensorname])
                val = _interp_pld_to_pld(pld, ds, val, indctd)
                ncvar['method'] = 'linear fill'
            else:
                val = gli[sensorname]
                #val = utils._zero_screen(val)
        #        val[val==0] = np.NaN
                val = convert(val)
                val = _interp_gli_to_pld(gli, ds, val, indctd)

            # make the attributes:
            ncvar[name].pop('coordinates', None)
            attrs = ncvar[name]
            attrs = utils.fill_required_attrs(attrs)
            ds[name] = (('time'), val, attrs)

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
    ds = utils.get_profiles(ds)

    ds = utils.get_derived_eos_raw(ds)

    ds = ds.assign_coords(longitude=ds.longitude)
    ds = ds.assign_coords(latitude=ds.latitude)
    ds = ds.assign_coords(depth=ds.depth)
    #ds = ds._get_distance_over_ground(ds)

    ds = utils.fill_metadata(ds, deployment['metadata'])
    start = ((ds['time'].values[0]).astype('timedelta64[s]') +
        np.datetime64('1970-01-01T00:00:00'))
    end = ((ds['time'].values[-1]).astype('timedelta64[s]')  +
        np.datetime64('1970-01-01T00:00:00'))

    ds.attrs['deployment_start'] = str(start)
    ds.attrs['deployment_end'] = str(end)

    try:
        os.mkdir('L1-timeseries')
    except:
        pass
    outname = 'L1-timeseries/' + ds.attrs['id'] +  '_L1.nc'
    _log.info('writing %s', outname)
    ds.to_netcdf(outname, 'w')

    return outname
