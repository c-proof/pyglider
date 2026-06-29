# -*- coding: utf-8 -*-
"""
SeaExplorer-specific processing routines.
"""

import datetime
import glob
import logging
import os
import warnings

import numpy as np
import polars as pl
import xarray as xr

import pyglider.utils as utils

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
    return outdir + fnout + 'parquet', filenum


def _needsupdating(ftype, fin, fout):
    if not os.path.isfile(fout):
        return True
    return os.path.getmtime(fin) >= os.path.getmtime(fout)


def _sort(ds):
    return ds.sortby('time')


def raw_to_rawnc(
    indir,
    outdir,
    deploymentyaml,
    incremental=True,
    min_samples_in_file=5,
    dropna_subset=None,
    dropna_thresh=1,
):
    """
    Convert seaexplorer text files to raw parquet pandas files.

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

    dropna_subset : list of strings, default None
        If more values than *dropna_thresh* of the variables listed here are
        empty (NaN), then drop this line of data.  Useful for raw payload files
        that are heavily oversampled.  Get the variable names from the raw text
        file.  See `pandas.DataFrame.dropna`.

    dropna_thresh : integer, default 1
        Number of variables listed in dropna_subset that can be empty before
        the line is dropped.


    Returns
    -------
    status : bool
        *True* success.

    Notes
    -----

    This process can be slow for many files.

    For the *dropna* functionality, list one variable for each of the sensors
    that is *not* over-sampled.  For instance, we had an AROD, GPCTD, and
    FLBBCD and the AROD was grossly oversampled, whereas the other two were not,
    but were not sampled synchronously.  In that case we chose:
    `dropna_subset=['GPCTD_TEMPERATURE', 'FLBBCD_CHL_COUNT']` to keep all
    rows where either of these were good, and dropped all other rows.

    """
    # Create out directory for netcdfs if it does not exist
    try:
        os.mkdir(outdir)
    except FileExistsError:
        pass

    for ftype in ['gli', 'pld1']:
        goodfiles = []
        badfiles = []
        for rawsub in ['raw', 'sub']:
            _log.info(f'Reading in raw files matching *{ftype}.{rawsub}*')
            d = indir + f'*.{ftype}.{rawsub}.*'

            files = glob.glob(d)
            _log.info(f'Nfiles found: {len(files)}')

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

            for ind, f in enumerate(files):
                # output name:
                fnout, filenum = _outputname(f, outdir)
                _log.info(f'{f} to {fnout}')
                if not incremental or _needsupdating(ftype, f, fnout):
                    _log.info(f'Doing: {f} to {fnout}')
                    # Try to read the file with polars. If the file is corrupted (rare), file read will fail and file
                    # is appended to badfiles
                    try:
                        out = pl.read_csv(f, separator=';')
                    except Exception as e:
                        _log.warning(f'Exception reading {f}: {e}')
                        _log.warning(f'Could not read {f}')
                        badfiles.append(f)
                        continue
                    # Parse the datetime from nav files (called Timestamp) and pld1 files (called PLD_REALTIMECLOCK)
                    if 'Timestamp' in out.columns:
                        out = out.with_columns(
                            pl.col('Timestamp').str.strptime(
                                pl.Datetime, format='%d/%m/%Y %H:%M:%S'
                            )
                        )
                        out = out.rename({'Timestamp': 'time'})
                    else:
                        out = out.with_columns(
                            pl.col('PLD_REALTIMECLOCK').str.strptime(
                                pl.Datetime, format='%d/%m/%Y %H:%M:%S.%3f'
                            )
                        )
                        out = out.rename({'PLD_REALTIMECLOCK': 'time'})
                    # strip leading and trailing spaces
                    out = out.with_columns(pl.col(pl.String).str.strip_chars())
                    # replace empty values with None
                    out = out.with_columns(pl.when(pl.col(pl.String).str.len_chars() == 0)
                                           .then(None)
                                           .otherwise(pl.col(pl.String))
                                           .name.keep())
                    for col_name in out.columns:
                        if 'time' not in col_name.lower():
                            out = out.with_columns(pl.col(col_name).cast(pl.Float64))
                    # remove leading and trailing spaces from column names
                    out = out.rename(str.strip)
                    # If AD2CP data present, convert timestamps to datetime
                    if 'AD2CP_TIME' in out.columns:
                        # Set datestamps with date 00000 to None
                        out = out.with_columns(
                            pl.col('AD2CP_TIME').str.strptime(
                                pl.Datetime, format='%m%d%y %H:%M:%S', strict=False
                            )
                        )

                    # subsetting for heavily oversampled raw data:
                    if rawsub == 'raw' and dropna_subset is not None:
                        # This check is the polars equivalent of pandas dropna. See docstring note on dropna
                        dropna_thresh = 1  # Drop rows with more than 1 null
                        null_exprs = [pl.col(c).is_null().cast(pl.Int64) for c in df.columns]

                        df = df.with_columns(
                            pl.sum(null_exprs).alias("nulls")
                        ).filter(
                            pl.col("nulls") <= dropna_thresh
                        ).drop("nulls")
                    if ftype == 'gli':
                        out = out.with_columns(
                            [(pl.col('NavState') * 0 + int(filenum)).alias('fnum')]
                        )
                        out.write_parquet(fnout)
                        goodfiles.append(f)
                    else:
                        if out.select('time').shape[0] > min_samples_in_file:
                            out.write_parquet(fnout)
                            goodfiles.append(f)
                        else:
                            _log.warning(
                                'Number of sensor data points'
                                'too small. Skipping parquet write'
                            )
                            badfiles.append(f)
    if len(badfiles) > 0:
        _log.warning('Some files could not be parsed:')
        for fn in badfiles:
            _log.warning('%s', fn)
    if not goodfiles:
        _log.warning(f'No valid unprocessed seaexplorer files found in' f'{indir}')
        return False
    _log.info('All raw files converted to parquet')
    return True


def drop_pre_1971_samples(df):
    """
    If dates greater than 1971, 1, 1 are observed, drop any dates before
    1971-01-01, from the datset and return it. This function removes 1970
    timestamps caused by a SeaExplorer rebooting during a mission. If all dates
    are < 1971-01-01, no action is taken

    Parameters:
        df: polars.DataFrame
            dataframe to check for pre-1971 dates
    Returns:
        df: polars.DataFrame
    """
    dt_1971 = datetime.datetime(1971, 1, 1)
    # If all dates before or after 1971-01-01, return the dataset
    pre_1971 = df.filter(pl.col('time') < dt_1971)
    if len(pre_1971) == len(df):
        return pre_1971
    post_1971 = df.filter(pl.col('time') > dt_1971)
    if len(post_1971) == len(df):
        return post_1971
    return df.filter(pl.col('time') > dt_1971)


def merge_parquet(indir, outdir, deploymentyaml, incremental=False, kind='raw'):
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

    deployment = utils._get_deployment(deploymentyaml)

    metadata = deployment['metadata']
    id = metadata['glider_name']
    outgli = outdir + '/' + id + '-rawgli.parquet'
    outpld = outdir + '/' + id + '-' + kind + 'pld.parquet'

    _log.info('Opening *.gli.sub.*.parquet multi-file dataset from %s', indir)
    files = sorted(glob.glob(indir + '/*.gli.sub.*.parquet'))
    if not files:
        _log.warning(f'No *gli*.parquet files found in {indir}')
        return False
    # lets figure out what columns to read for situations where the number of columns changes:
    gli0 = pl.read_parquet(files[0])
    gli1 = pl.read_parquet(files[-1])
    columns = list(set(gli0.columns) | set(gli1.columns))
    columns = set([c for c in columns if len(c) > 0])
    dfs = []
    for f in files:
        df = pl.read_parquet(f)
        missing = columns - set(df.columns)
        for col in missing:
            df = df.with_columns(pl.lit(None).cast(pl.Float64).alias(col))

        dfs.append(df.select(sorted(columns)))
    gli = pl.concat(dfs, rechunk=True)

    #gli = pl.read_parquet(indir + '/*.gli.sub.*.parquet', columns=columns,
    #                      missing_columns='insert')
    gli = drop_pre_1971_samples(gli)
    gli.write_parquet(outgli)
    _log.info(f'Done writing {outgli}')

    _log.info('Opening *.pld.sub.*.parquet multi-file dataset')
    files = sorted(glob.glob(indir + '/*.pld1.' + kind + '.*.parquet'))
    if not files:
        _log.warning(f'No *{kind}*.parquet files found in {indir}')
        return False
    gli0 = pl.read_parquet(files[0])
    gli1 = pl.read_parquet(files[-1])
    columns = list(set(gli0.columns) | set(gli1.columns))
    columns = set([c for c in columns if len(c) > 0])
    dfs = []
    for f in files:
        df = pl.read_parquet(f)
        missing = columns - set(df.columns)
        for col in missing:
            df = df.with_columns(pl.lit(None).cast(pl.Float64).alias(col))

        dfs.append(df.select(sorted(columns)))
    pld = pl.concat(dfs, rechunk=True)

    #    pld = pl.read_parquet(indir + '/*.pld1.' + kind + '.*.parquet',
    #                          missing_columns='insert', columns=columns)
    pld = drop_pre_1971_samples(pld)
    pld.write_parquet(outpld)

    _log.info(f'Done writing {outpld}')
    _log.info('Done merge_rawnc')
    return True


def _interp_gli_to_pld(gli, ds, val, indctd):
    gli_ind = ~np.isnan(val)
    # switch for if we are comparing two polars dataframes or a polars dataframe and a xarray dataset
    if type(ds) is pl.DataFrame:
        valout = np.interp(ds['time'], gli.filter(gli_ind)['time'], val[gli_ind])
    else:
        valout = np.interp(
            ds['time'].astype(int),
            np.array(
                gli.filter(gli_ind)['time'].to_numpy().astype('datetime64[ns]')
            ).astype(int),
            val[gli_ind],
        )
    return valout


def _interp_pld_to_pld(pld, ds, val, indctd):
    pld_ind = np.where(~np.isnan(val))[0]
    if len(pld_ind) != len(indctd):
        val = np.interp(ds['time'], pld['time'][pld_ind], val[pld_ind])
    else:
        val = val[indctd]
    return val


def _remove_fill_values(df, fill_value=9999):
    """
    For input dataframe df, this function converts all Float values equaling fill_values to null. Columns of other
    datatypes are not affected.
    """
    df = df.with_columns(
        pl.when(pl.col(pl.Float64) == fill_value)
        .then(None)
        .otherwise(pl.col(pl.Float64))
        .name.keep()
    )
    return df


def _forward_fill(gli, todo='Lat'):
    """Forward-fill the specified column (todo) to propagate the last good value at each row."""
    gli = gli.with_columns([
        pl.col(todo).fill_null(strategy="forward").alias("temp_fill")
    ])
    gli = gli.with_columns([
        pl.when(
            (pl.col(todo) == pl.col("temp_fill").shift(1)) & pl.col(todo).is_not_null()
        ).then(np.nan).otherwise(pl.col(todo)).alias(todo)
    ])
    gli = gli.drop("temp_fill")
    return gli


def _drop_if(gli, todo='Lat', condit='DeadReckoning', value=1):
    """Drop Lat if DeadReckoning is 1"""
    gli = gli.with_columns([
        pl.when(pl.col(condit) == value).then(np.nan).otherwise(pl.col(todo)).alias(todo)
    ])
    return gli



def raw_to_timeseries(
    indir,
    outdir,
    deploymentyaml,
    kind='raw',
    profile_filt_time=100,
    profile_min_time=300,
    maxgap=10,
    interpolate=False,
    start_time=None,
    fnamesuffix='',
    deadreckon=False,
    replace_attrs=None
):
    """
    Convert raw seaexplorer data to a timeseries netcdf file.

    Parameters
    ----------
    indir : str
        Directory with the raw files are kept.

    outdir : str
        Directory to write the matching ``*.nc`` files.

    deploymentyaml : str
        YAML text file with deployment information for this glider.

    kind : 'raw' or 'sub'
        The type of data to process.  'raw' is the full resolution data, 'sub'
        is the sub-sampled data.  The default is 'raw'.  Note that realtime data is
        typically sub-sampled.

    profile_filt_time : float
        Time in seconds to use for filtering the profiles.  Default is 100.

    profile_min_time : float
        Minimum time in seconds for a profile to be considered a valid profile.
        Default is 300.

    maxgap : float
        Maximum gap in seconds to interpolate over.  Default is 10.

    interpolate : bool
        If *True*, interpolate the data to fill in gaps.  Default is False.

    start_time : str or None
        Drop data if before this date - sometimes there are bad times. Default is *None*

    fnamesuffix : str
        Suffix to add to the output file name.  Default is ''.

    deadreckon : bool
        If *True* use the dead reckoning latitude and longitude data from the glider.  Default
        is *False*, and latitude and longitude are linearly interpolated between surface fixes.
        *False* is the default, and recommended to avoid a-physical underwater jumps.

    replace_attrs : dict or None
        replace global attributes in the metadata after reading the metadata
        file in.  Helpful when processing runs with only a couple things that
        change.


    Returns
    -------
    outname : str
        Name of the output netcdf file.

    """

    deployment = utils._get_deployment(deploymentyaml)
    if replace_attrs:
        for att in replace_attrs:
            deployment['metadata'][att] = replace_attrs[att]

    metadata = deployment['metadata']
    ncvar = deployment['netcdf_variables']
    device_data = deployment['glider_devices']
    varnames = utils._get_varnames(deployment)
    time_name = varnames.get('time', 'time')
    lat_name = varnames.get('latitude', 'latitude')
    lon_name = varnames.get('longitude', 'longitude')
    pressure_name = varnames.get('pressure', 'pressure')
    depth_name = varnames.get('depth', 'depth')
    id = metadata['glider_name']
    _log.info(f'Opening combined nav file {indir}/{id}-rawgli.nc')
    gli = pl.read_parquet(f'{indir}/{id}-rawgli.parquet')
    _log.info(f'Opening combined payload file {indir}/{id}-{kind}pld.parquet')
    sensor = pl.read_parquet(f'{indir}/{id}-{kind}pld.parquet')
    sensor = _remove_fill_values(sensor)

    # don't use lat/lon if deadreckoned:
    if not deadreckon:
        lat_source = ncvar.get(lat_name, {}).get('source', '')
        lon_source = ncvar.get(lon_name, {}).get('source', '')
        if lat_source not in ('Lat', 'NAV_LATITUDE'):
            warnings.warn("For deadreckon=False, it is suggested to use 'Lat' or 'NAV_LATITUDE' as the source for latitude.")
        if lon_source not in ('Lon', 'NAV_LONGITUDE'):
            warnings.warn("For deadreckon=False, it is suggested to use 'Lon' or 'NAV_LONGITUDE' as the source for longitude.")
        if 'DeadReckoning' in gli.columns:
            _log.info('Not using deadreckoning; glider has DeadReckoning column')
            gli = _drop_if(gli, todo='Lat', condit='DeadReckoning', value=1)
            gli = _drop_if(gli, todo='Lon', condit='DeadReckoning', value=1)
        else:
            _log.info('Not using deadreckoning; glider does not have DeadReckoning column')
            gli = _drop_if(gli, todo='Lat', condit='NavState', value=116)
            gli = _drop_if(gli, todo='Lon', condit='NavState', value=116)
        # drop a lat/lon if it is not unique.  Happens when there
        # are stale fixes.
        gli = _forward_fill(gli, todo='Lat')
        gli = _forward_fill(gli, todo='Lon')

    # build a new data set based on info in `deploymentyaml.`
    # We will use ctd as the interpolant
    ds = xr.Dataset()
    attr = {}
    for atts in ncvar[time_name].keys():
        if atts != 'coordinates':
            attr[atts] = ncvar[time_name][atts]

    # If present, use the timebase specified in ncvar: timebase in the
    # deployment yaml.
    if 'timebase' not in ncvar:
        raise ValueError(
            'Must specify timebase:source in netcdf_variables section of deployment yaml'
        )
    if ncvar['timebase']['source'] not in sensor.columns:
        raise ValueError(
            f"timebase source: {ncvar['timebase']['source']} not found in pld1 columns"
        )
    vals = sensor.select([ncvar['timebase']['source']]).to_numpy()[:, 0]
    indctd = np.where(~np.isnan(vals))[0]
    ds['time'] = (
        ('time'),
        sensor.select('time').to_numpy()[indctd, 0].astype('datetime64[ns]'),
        attr,
    )
    thenames = list(ncvar.keys())
    # Check yaml to see if interpolate has been set to True
    if 'interpolate' in thenames:
        if ncvar['interpolate']:
            interpolate = True
    for i in [time_name, 'timebase', 'keep_variables', 'interpolate']:
        if i in thenames:
            thenames.remove(i)
    for name in thenames:
        _log.info('interpolating ' + name)
        if ('method' not in ncvar[name].keys()
                and 'processing_method' not in ncvar[name].keys()):
            if 'source' not in ncvar[name]:
                # no source and no processing_method — skip (e.g. QC placeholders)
                continue
            # variables that are in the data set or can be interpolated from it
            if 'conversion' in ncvar[name].keys():
                convert = getattr(utils, ncvar[name]['conversion'])
            else:
                convert = utils._passthrough
            sensorname = ncvar[name]['source']
            if sensorname in list(sensor.columns):
                _log.debug('sensorname %s', sensorname)
                val = convert(sensor.select(sensorname).to_numpy()[:, 0])
                if 'coarsen' in ncvar[name]:
                    # coarsen oxygen data as originally perscribed
                    coarsen_time = ncvar[name]['coarsen']
                    # create a boolean mask of coarsened timesteps. Use this mask to create an array of samples to keep
                    coarse_ints = np.arange(
                        0, len(sensor) / coarsen_time, 1 / coarsen_time
                    ).astype(int)
                    sensor_sub = sensor.with_columns(
                        pl.lit(coarse_ints).alias('coarse_ints')
                    )
                    # Subsample the variable data keeping only the samples from the coarsened timeseries
                    sensor_sub_grouped = (
                        sensor_sub.with_columns(pl.col('time').to_physical())
                        .group_by(pl.col('coarse_ints'), maintain_order=True)
                        .mean()
                        .with_columns(pl.col('time').cast(pl.Datetime('ms')))[:-1, :]
                    )
                    val2 = sensor_sub_grouped.select(sensorname).to_numpy()[:, 0]
                    val = _interp_gli_to_pld(sensor_sub_grouped, sensor, val2, indctd)
                if interpolate and not np.isnan(val).all():
                    time_original = sensor.select('time').to_numpy()[:, 0]
                    time_var = time_original[np.where(~np.isnan(val))[0]]
                    var_non_nan = val[np.where(~np.isnan(val))[0]]
                    time_timebase = sensor.select('time').to_numpy()[indctd, 0]
                    if val.dtype == '<M8[us]':
                        # for datetime, must convert to numerical, interpolate, then convert back
                        us_since_1970 = (
                            var_non_nan - np.datetime64('1970-01-01')
                        ).astype(int)
                        val_int = np.interp(
                            time_timebase.astype(float),
                            time_var.astype(float),
                            us_since_1970,
                        )
                        val_us = val_int.astype('timedelta64[us]')
                        val = np.datetime64('1970-01-01') + val_us
                    else:
                        val = np.interp(
                            time_timebase.astype(float),
                            time_var.astype(float),
                            var_non_nan,
                        )

                    # interpolate only over those gaps that are smaller than 'maxgap'
                    # apparently maxgap is to be in somethng like seconds, and this data is in ms.  Certainly
                    # the default of 0.3 s was not intended.  Changing default to 10 s:
                    tg_ind = utils.find_gaps(
                        time_var.astype(float),
                        time_timebase.astype(float),
                        maxgap * 1000,
                    )
                    val[tg_ind] = np.nan
                else:
                    val = val[indctd]

                ncvar['method'] = 'linear fill'
            else:
                val = gli.select(sensorname).to_numpy()[:, 0]
                val = convert(val)
                # Values from the glider netcdf must be interpolated to match
                # the sensor netcdf
                val = _interp_gli_to_pld(gli, ds, val, indctd)

            # make the attributes:
            ncvar[name].pop('coordinates', None)
            attrs = ncvar[name]
            attrs = utils.fill_required_attrs(attrs)
            ds[name] = (('time'), val, attrs)

    # fix lon and lat to be linearly interpolated between fixes
    lat_name = utils._resolve_role(ds, varnames, 'latitude')
    lon_name = utils._resolve_role(ds, varnames, 'longitude')
    good = (
        np.where(
            np.abs(np.diff(ds[lon_name])) + np.abs(np.diff(ds[lat_name])) > 0
        )[0]
        + 1
    )
    ds[lon_name].values = np.interp(ds.time, ds.time[good], ds[lon_name][good])
    ds[lat_name].values = np.interp(ds.time, ds.time[good], ds[lat_name][good])

    # keep only timestamps with data from one of a set of variables
    if 'keep_variables' in ncvar:
        keeps = np.empty(len(ds[lon_name]))
        keeps[:] = np.nan
        keeper_vars = ncvar['keep_variables']
        for keep_var in keeper_vars:
            keeps[~np.isnan(ds[keep_var].values)] = 1
        ds = ds.where(~np.isnan(keeps))
        ds = ds.dropna(dim='time', how='all')

    # drop dates before start_time
    if start_time is not None:
        ds = ds.where(ds.time >= np.datetime64(start_time), drop=True)

    # Derived variables — depth, profiles, distance, thermodynamics.
    # For OG 1.0 YAMLs these are specified via processing_method; for IOOS GDAC
    # YAMLs we fall back to the legacy utility functions.
    has_dog_method = any(
        isinstance(a, dict) and 'processing_method' in a
        and 'distance_over_ground' in a['processing_method']
        for a in ncvar.values()
    )
    has_thermo_methods = any(
        isinstance(a, dict) and 'processing_method' in a
        and any(m in utils._THERMO_METHODS for m in a['processing_method'])
        for a in ncvar.values()
    )

    if pressure_name in ds:
        ds = utils.get_glider_depth(ds, varnames=varnames)
    if not has_dog_method:
        ds = utils.get_distance_over_ground(ds, varnames=varnames)
    ds = utils.get_profiles_new(
        ds, filt_time=profile_filt_time, profile_min_time=profile_min_time,
        varnames=varnames,
    )
    if not has_thermo_methods:
        cond_name = varnames.get('conductivity', 'conductivity')
        temp_name = varnames.get('temperature', 'temperature')
        if all(n in ds for n in (temp_name, cond_name, pressure_name)):
            ds = utils.get_derived_eos_raw(ds, varnames=varnames)

    # Dispatch processing_method derived variables (OG 1.0: PSAL, SIGMA0,
    # DEPTH, DISTANCE_OVER_GROUND, PROFILE_NUMBER, …).
    ds = utils._dispatch_processing_methods(ds, ncvar)

    # somehow this comes out unsorted:
    ds = ds.sortby(ds.time)
    # Drop duplicate timestamps and check how many are removed this way
    len_before_drop = len(ds.time)
    if hasattr(ds, 'drop_duplicates'):
        ds = ds.drop_duplicates(dim='time')
    else:
        time_diffs = (ds.time.astype(int).diff(dim='time') > 1e-6).values
        time_diffs_list = list(time_diffs)
        time_diffs_list.append(True)
        good = np.array(time_diffs_list)
        ds = ds.isel(time=good)
    len_after_drop = len(ds.time)
    proportion_kept = len_after_drop / len_before_drop
    loss_str = (
        f'{100 * (1 - proportion_kept)} % samples removed by timestamp deduplication.'
    )
    if proportion_kept < 0.5:
        raise ValueError(f'{loss_str} Check input data for duplicate timestamps')
    elif proportion_kept < 0.999:
        _log.warning(loss_str)
    else:
        _log.info(loss_str)

    # Correct oxygen if present (IOOS GDAC format only — looks for variable
    # named 'oxygen_concentration' with a 'correct_oxygen' sub-key):
    if 'oxygen_concentration' in ncvar.keys():
        if 'correct_oxygen' in ncvar['oxygen_concentration'].keys():
            ds = utils.oxygen_concentration_correction(ds, ncvar)
        else:
            _log.warning(
                'correct_oxygen not found in oxygen yaml. No correction applied'
            )

    if lon_name in ds:
        ds = ds.assign_coords({lon_name: ds[lon_name]})
    if lat_name in ds:
        ds = ds.assign_coords({lat_name: ds[lat_name]})
    if depth_name in ds:
        ds = ds.assign_coords({depth_name: ds[depth_name]})

    ds = utils.fill_metadata(ds, deployment['metadata'], device_data,
                             varnames=varnames)

    # GPS fix variables: compute once, assign to all variables that specify
    # processing_method: gps_fixes_from_nav.  Source columns and per-variable
    # attributes come entirely from the YAML — no hard-coded trigger here.
    gps_var_names = [
        name for name, attrs in ncvar.items()
        if isinstance(attrs, dict)
        and isinstance(attrs.get('processing_method'), dict)
        and 'gps_fixes_from_nav' in attrs['processing_method']
    ]
    if gps_var_names:
        try:
            # Use source columns from the first entry (they must all agree).
            first_inputs = ncvar[gps_var_names[0]]['processing_method']['gps_fixes_from_nav'] or {}
            lat_col = first_inputs.get('lat_source', 'Lat')
            lon_col = first_inputs.get('lon_source', 'Lon')

            gli_raw = pl.read_parquet(f'{indir}/{id}-rawgli.parquet')
            if 'DeadReckoning' in gli_raw.columns:
                gli_gps = gli_raw.filter(pl.col('DeadReckoning') == 0)
            else:
                gli_gps = gli_raw.filter(pl.col('NavState') != 116)
            gli_gps = gli_gps.drop_nulls(subset=[lat_col, lon_col])
            lat_gps = utils.nmea2deg(gli_gps[lat_col].to_numpy())
            lon_gps = utils.nmea2deg(gli_gps[lon_col].to_numpy())
            t_gps = gli_gps['time'].to_numpy().astype('datetime64[ns]')
            valid = np.isfinite(lat_gps) & np.isfinite(lon_gps) & (lat_gps != 0) & (lon_gps != 0)
            lat_gps, lon_gps, t_gps = lat_gps[valid], lon_gps[valid], t_gps[valid]

            # Map GPS fixes onto the sensor time grid (NaN elsewhere).
            n = len(ds['time'])
            ds_times_ns = ds['time'].values.astype(np.int64)
            gps_times_ns = t_gps.astype(np.int64)
            lat_out = np.full(n, np.nan)
            lon_out = np.full(n, np.nan)
            t_out = np.full(n, np.nan)
            for i in range(len(gps_times_ns)):
                idx = np.searchsorted(ds_times_ns, gps_times_ns[i])
                if idx >= n:
                    idx = n - 1
                elif idx > 0 and (abs(ds_times_ns[idx - 1] - gps_times_ns[i]) <
                                  abs(ds_times_ns[idx] - gps_times_ns[i])):
                    idx -= 1
                lat_out[idx] = lat_gps[i]
                lon_out[idx] = lon_gps[i]
                t_out[idx] = gps_times_ns[i] / 1e9  # seconds since 1970-01-01

            role_data = {'latitude': lat_out, 'longitude': lon_out, 'time': t_out}

            for varname in gps_var_names:
                inputs = ncvar[varname]['processing_method']['gps_fixes_from_nav'] or {}
                role = inputs.get('role')
                if role not in role_data:
                    _log.warning('gps_fixes_from_nav: unknown role %r for %s; skipping', role, varname)
                    continue
                var_attrs = {k: v for k, v in ncvar[varname].items()
                             if k not in ('processing_method', 'processing_role',
                                          'average_method', 'source', 'coordinates')}
                var_attrs = utils.fill_required_attrs(var_attrs)
                ds[varname] = (('time',), role_data[role], var_attrs)

            n_gps = int(np.sum(np.isfinite(lat_out)))
            _log.info('Added %d GPS fixes via gps_fixes_from_nav for %s',
                      n_gps, gps_var_names)
        except Exception:
            _log.warning('Could not extract GPS fix variables from gli data',
                         exc_info=True)

    start = ds['time'].values[0]
    end = ds['time'].values[-1]

    ds.attrs['deployment_start'] = str(start)
    ds.attrs['deployment_end'] = str(end)

    try:
        os.mkdir(outdir)
    except:
        pass
    id0 = ds.attrs['deployment_name']
    outname = outdir + id0 + fnamesuffix + '.nc'
    _log.info('writing %s', outname)
    if 'units' in ds.time.attrs.keys():
        ds.time.attrs.pop('units')
    if 'calendar' in ds.time.attrs.keys():
        ds.time.attrs.pop('calendar')
    if 'ad2cp_time' in list(ds):
        if 'units' in ds.ad2cp_time.attrs.keys():
            ds.ad2cp_time.attrs.pop('units')
    ds = utils.make_scalar_variables(ds, deployment)
    ds = utils.make_sensor_variables(ds, deployment)
    utils._save_dataset(
        ds,
        outname,
        deployment,
        mode='w',
        encoding={
            'time': {'units': 'seconds since 1970-01-01T00:00:00Z', 'dtype': 'float64'}
        },
    )
    return outname


def parse_surfacing_logfile(fn):
    """
    Parse the surfacing log file from a seaexplorer deployment.  This is not
    used in the processing, but is useful to track the glider behaviour if nav
    data not being sent back.

    Parameters
    ----------
    fn : str
        Path to the log file.

    Returns
    -------
    df : polars.DataFrame
        Dataframe with the log data.

    """

    """
        "18.05.2026 01:07:07.296";TRACE;"module-IrisCom";"Glider";"SEA035";"done (a): Waiting for GO"
        "18.05.2026 01:07:07.499";TRACE;"module-IrisCom";"Glider";"SEA035";"$SEAMRS,SEA035,143,273,0,275,180526,010709,4825.877,-12522.708*1106;"
        "18.05.2026 01:07:07.499";TRACE;"module-IrisCom";"Glider";"SEA035";"$SEANAV,1,41,1020,-1020,1020,-1020,1000,2,10,10*5a46;"
        "18.05.2026 01:07:07.702";TRACE;"module-IrisCom";"Glider";"SEA035";"$SEADST,116,279,16,-48,-7,-0.3,76078,9.6,275,100,0,500*fdfa;"
        "18.05.2026 01:07:10.718";TRACE;"module-IrisCom";"Glider";"SEA035";"$SEALOG,1,0,148.4,75897,78006,8.3,9.6,27.1,27.7,-0.18,0.18,-19.4,19.4,44.8,30.3,-4
        78.1,260.1*db10;"
        "18.05.2026 01:07:10.718";TRACE;"module-IrisCom";"Glider";"SEA035";"$SEADEV,1,1,1,1,1*2c63;"
        "18.05.2026 01:07:10.718";TRACE;"module-IrisCom";"Glider";"SEA035";"$SEADRK,0,4825.877,-12522.708,0.00,3142,0.22,113*5d1b;"
        "18.05.2026 01:07:27.374";TRACE;"module-IrisCom";"Glider";"SEA035";"done (a): Waiting for GO"
        "18.05.2026 01:07:30.077";TRACE;"module-IrisCom";"Glider";"SEA035";"$SEAMRS,SEA035,143,273,0,276,180526,010729,4825.873,-12522.705*5ad5;"
        "18.05.2026 01:07:30.077";TRACE;"module-IrisCom";"Glider";"SEA035";"$SEANAV,1,41,1020,-1020,1020,-1020,1000,2,10,10*5a46;"

        Translated as
        $SEAMRS,<ID>,<MSN>,<CYCLE>,<STATUS>,<BATTERIES>,<DATE>,<TIME>,<LAT>,<LON>*<CKS>;
        $SEANAV,<MODE>,<HEADING>,<PU>,<PD>,<BU>,<BD>,<ZB>,<ZT>,<AL>,<RATE>*<CKS>;
        $SEADST,<NAVSTATE>,<HEADING>,<DECLINATION>,<PITCH>,<ROLL>,<DEPTH>,<VACUUM>,<TEMP>,<VOLTAGE>,<LPOS>,<APOS>,<BPOS>*<CKS>;
        $SEALOG,<CALLS>,<CALLS_FAILED>,<ZMAX>,<VAC_MIN>,<VAC_MAX>,<TEMP_MIN>,<TEMP_MAX>,<VOLT_MIN>,<VOLT_MAX>,<VEL_D>,<VEL_A>,<PITCH_D>,<PITCH_A>,<LIN_D>,<LIN_A>,<BAL_D>,<BAL_A>*<CKS>;
        $SEADEV,<RADIO>,<GPS>,<IRIDIUM>*<CKS>
    """

    # we need to go through this file and extract from each line a time series of
    # time, lat, lon, battery and heading commanded (from SEANAV):
    # each SEAMRS should have a correspondiong SEANAV so use the SEAMRS as the index and
    # pull out the time, lat, lon, battery from it, and then pull out the heading from
    # the corresponding SEANAV.  We can also pull out the NAVSTATE from SEADST if we want
    # to know when the glider is dead reckoning.

    df = pl.read_csv(
        fn,
        separator=';',
        has_header=False,
        new_columns=[
            'timestamp',
            'log_level',
            'module',
            'device_type',
            'device_id',
            'message',
        ],
    )
    df = df.with_columns(
        pl.col('timestamp').str.strptime(
            pl.Datetime,
            format='%d.%m.%Y %H:%M:%S.%3f',
            strict=False,
        )
    )
    df = df.rename({'timestamp': 'time'})
    df = df.filter(pl.col('time').is_not_null())

    # Parse only Iris/SEA NMEA-like payloads in the message column.
    parsed = (
        df.with_row_count('row_id')
        .with_columns(
            pl.col('message')
            .str.strip_chars()
            .str.strip_prefix('"')
            .str.strip_suffix('"')
            .str.strip_suffix(';')
            .str.strip_prefix('$')
            .alias('message_clean')
        )
        .with_columns(
            pl.col('message_clean').str.extract(r'\*([0-9A-Fa-f]+)$', 1).alias('checksum'),
            pl.col('message_clean')
            .str.replace(r'\*[0-9A-Fa-f]+$', '')
            .alias('payload')
        )
        .filter(pl.col('payload').str.starts_with('SEA'))
        .with_columns(pl.col('payload').str.split(',').alias('parts'))
        .with_columns(pl.col('parts').list.get(0).alias('record_type'))
    )

    seamrs = parsed.filter(
        (pl.col('record_type') == 'SEAMRS') & (pl.col('parts').list.len() >= 10)
    ).select(
        'row_id',
        'time',
        pl.col('checksum').alias('seamrs_checksum'),
        pl.col('parts').list.get(1).alias('id'),
        pl.col('parts').list.get(2).cast(pl.Int64, strict=False).alias('msn'),
        pl.col('parts').list.get(3).cast(pl.Int64, strict=False).alias('cycle'),
        pl.col('parts').list.get(4).cast(pl.Int64, strict=False).alias('status'),
        (
            pl.col('parts').list.get(5).cast(pl.Float64, strict=False) / 10.0
        ).alias('battery'),
        pl.col('parts').list.get(8).cast(pl.Float64, strict=False).alias('lat'),
        pl.col('parts').list.get(9).cast(pl.Float64, strict=False).alias('lon'),
    )



    if seamrs.is_empty():
        return seamrs

    seanav = parsed.filter(
        (pl.col('record_type') == 'SEANAV') & (pl.col('parts').list.len() >= 5)
    ).select(
        'row_id',
        pl.col('checksum').alias('seanav_checksum'),
        pl.col('parts').list.get(1, null_on_oob=True).cast(pl.Int64, strict=False).alias('mode'),
        pl.col('parts').list.get(2, null_on_oob=True).cast(pl.Float64, strict=False).alias('heading_cmd'),
        pl.col('parts').list.get(3, null_on_oob=True).cast(pl.Float64, strict=False).alias('pu'),
        pl.col('parts').list.get(4, null_on_oob=True).cast(pl.Float64, strict=False).alias('pd'),
        pl.col('parts').list.get(5, null_on_oob=True).cast(pl.Float64, strict=False).alias('bu'),
        pl.col('parts').list.get(6, null_on_oob=True).cast(pl.Float64, strict=False).alias('bd'),
        pl.col('parts').list.get(8, null_on_oob=True).cast(pl.Float64, strict=False).alias('zt'),
        pl.col('parts').list.get(9, null_on_oob=True).cast(pl.Float64, strict=False).alias('al'),
        pl.col('parts').list.get(10, null_on_oob=True).cast(pl.Float64, strict=False).alias('rate'),
    )

    seadst = parsed.filter(
        (pl.col('record_type') == 'SEADST') & (pl.col('parts').list.len() >= 2)
    ).select(
        'row_id',
        pl.col('checksum').alias('seadst_checksum'),
        pl.col('parts').list.get(1).cast(pl.Int64, strict=False).alias('navstate'),
    )

    out = seamrs.sort('row_id')
    if not seanav.is_empty():
        out = out.join_asof(
            seanav.sort('row_id'),
            on='row_id',
            strategy='forward',
        )
    if not seadst.is_empty():
        out = out.join_asof(
            seadst.sort('row_id'),
            on='row_id',
            strategy='forward',
        )

    out = out.drop('row_id').sort('time')

    # Lat and Lon are NMEA lat and lon so need to be converted to decimal degrees:
    out = out.with_columns(
        pl.col('lat').map_elements(lambda x: utils.nmea2deg(x) if x is not None else None, return_dtype=pl.Float64).alias('lat'),
        pl.col('lon').map_elements(lambda x: utils.nmea2deg(x) if x is not None else None, return_dtype=pl.Float64).alias('lon'),
    )

    print(out[['time', 'lat', 'lon']])
    return out

def add_latlon_to_gridfiles(outname, df):
    """
    Add lat and lon to the grid files.  This is needed for the glider toolbox
    gridding to work.

    Parameters
    ----------
    outname : str
        Path to the output netcdf file.

    df : polars.DataFrame
        Dataframe with the lat and lon data.  Must have columns 'time', 'lat',
        and 'lon'.
    """
    with xr.open_dataset(outname, mode='r') as ds:
        ds.load()
        # interpolate lat and lon to the grid file time steps
        # Filter to only good (non-NaN) lat/lon values
        df_good = df.filter((pl.col('lat').is_not_null()) & (pl.col('lon').is_not_null()))

        # Convert both time arrays to same units (nanoseconds) for interpolation
        ds_time_ns = ds['time'].values.astype('datetime64[ns]').astype(float)
        df_time_ns = df_good['time'].to_numpy().astype('datetime64[ns]').astype(float)

        lat_interp = np.interp(
            ds_time_ns,
            df_time_ns,
            df_good['lat'].to_numpy(),
        )
        lon_interp = np.interp(
            ds_time_ns,
            df_time_ns,
            df_good['lon'].to_numpy(),
        )

        # Replace latitude and longitude only where they match the previous value
        # (indicating forward-filled bad values)
        lat_orig = ds['latitude'].values.copy()
        lon_orig = ds['longitude'].values.copy()

        # Find where current value matches previous value within tolerance
        # (forward-filled bad values can differ by tiny roundoff amounts).
        tol = 1e-6
        lat_repeated = np.concatenate(
            ([False], np.isclose(lat_orig[1:], lat_orig[:-1], atol=tol, rtol=0.0))
        )
        lon_repeated = np.concatenate(
            ([False], np.isclose(lon_orig[1:], lon_orig[:-1], atol=tol, rtol=0.0))
        )

        # Create new arrays, replacing only the repeated values
        lat_new = lat_orig.copy()
        lon_new = lon_orig.copy()
        lat_new[lat_repeated] = lat_interp[lat_repeated]
        lon_new[lon_repeated] = lon_interp[lon_repeated]

        ds['latitude'].values = lat_new
        ds['longitude'].values = lon_new

    ds.to_netcdf(outname, mode='w')



# aliases:
raw_to_L1timeseries = raw_to_L0timeseries = raw_to_timeseries
merge_rawnc = merge_parquet

__all__ = ['raw_to_rawnc', 'merge_parquet', 'raw_to_timeseries']
