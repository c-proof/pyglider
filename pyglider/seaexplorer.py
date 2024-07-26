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
import datetime
import polars as pl

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
    return (os.path.getmtime(fin) >= os.path.getmtime(fout))


def _sort(ds):
    return ds.sortby('time')


def raw_to_rawnc(indir, outdir, deploymentyaml, incremental=True,
                 min_samples_in_file=5, dropna_subset=None, dropna_thresh=1):
    """
    Convert seaexplorer text files to raw parquet files.

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
                    if "Timestamp" in out.columns:
                        out = out.with_columns(
                            pl.col("Timestamp").str.strptime(pl.Datetime, format="%d/%m/%Y %H:%M:%S"))
                        out = out.rename({"Timestamp": "time"})
                    else:
                        out = out.with_columns(
                            pl.col("PLD_REALTIMECLOCK").str.strptime(pl.Datetime, format="%d/%m/%Y %H:%M:%S.%3f"))
                        out = out.rename({"PLD_REALTIMECLOCK": "time"})
                    for col_name in out.columns:
                        if "time" not in col_name.lower():
                            out = out.with_columns(pl.col(col_name).cast(pl.Float64))
                    # If AD2CP data present, convert timestamps to datetime
                    if 'AD2CP_TIME' in out.columns:
                        # Set datestamps with date 00000 to None
                        out = out.with_columns(
                            pl.col('AD2CP_TIME').str.strptime(pl.Datetime, format="%m%d%y %H:%M:%S", strict=False))

                    # subsetting for heavily oversampled raw data:
                    if rawsub == 'raw' and dropna_subset is not None:
                        # This check is the polars equivalent of pandas dropna. See docstring note on dropna
                        out = out.with_columns(out.select(pl.col(dropna_subset).is_null().cast(pl.Int64))
                                              .sum(axis=1).alias("null_count")).filter(
                            pl.col("null_count") <= dropna_thresh) \
                            .drop("null_count")

                    if ftype == 'gli':
                        out = out.with_columns([(pl.col("NavState") * 0 + int(filenum)).alias("fnum")])
                        out.write_parquet(fnout)
                        goodfiles.append(f)
                    else:
                        if out.select("time").shape[0] > min_samples_in_file:
                            out.write_parquet(fnout)
                            goodfiles.append(f)
                        else:
                            _log.warning('Number of sensor data points'
                                         'too small. Skipping parquet write')
                            badfiles.append(f)
    if len(badfiles) > 0:
        _log.warning('Some files could not be parsed:')
        for fn in badfiles:
            _log.warning('%s', fn)
    if not goodfiles:
        _log.warning(f'No valid unprocessed seaexplorer files found in'f'{indir}')
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
    pre_1971 = df.filter(pl.col("time") < dt_1971)
    if len(pre_1971) == len(df):
        return pre_1971
    post_1971 = df.filter(pl.col("time") > dt_1971)
    if len(post_1971) == len(df):
        return post_1971
    return df.filter(pl.col("time") > dt_1971)


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
    gli = pl.read_parquet(indir + '/*.gli.sub.*.parquet')
    gli = drop_pre_1971_samples(gli)
    gli.write_parquet(outgli)
    _log.info(f'Done writing {outgli}')

    _log.info('Opening *.pld.sub.*.parquet multi-file dataset')
    files = sorted(glob.glob(indir + '/*.pld1.' + kind + '.*.parquet'))
    if not files:
        _log.warning(f'No *{kind}*.parquet files found in {indir}')
        return False
    pld = pl.read_parquet(indir + '/*.pld1.' + kind + '.*.parquet')
    pld = drop_pre_1971_samples(pld)
    pld.write_parquet(outpld)

    _log.info(f'Done writing {outpld}')
    _log.info('Done merge_rawnc')
    return True


def _interp_gli_to_pld(gli, ds, val, indctd):
    gli_ind = ~np.isnan(val)
    # switch for if we are comparing two polars dataframes or a polars dataframe and a xarray dataset
    if type(ds) is pl.DataFrame:
        valout = np.interp(ds["time"],
                           gli.filter(gli_ind)["time"],
                           val[gli_ind])
    else:
        valout = np.interp(ds['time'].astype(int),
                           np.array(gli.filter(gli_ind)["time"].to_numpy().astype('datetime64[ns]')).astype(int),
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


def raw_to_timeseries(indir, outdir, deploymentyaml, kind='raw',
                      profile_filt_time=100, profile_min_time=300,
                      maxgap=10, interpolate=False, fnamesuffix=''):
    """
    A little different than above, for the 4-file version of the data set.
    """

    deployment = utils._get_deployment(deploymentyaml)

    metadata = deployment['metadata']
    ncvar = deployment['netcdf_variables']
    device_data = deployment['glider_devices']
    id = metadata['glider_name']
    _log.info(f'Opening combined nav file {indir}/{id}-rawgli.nc')
    gli = pl.read_parquet(f'{indir}/{id}-rawgli.parquet')
    _log.info(f'Opening combined payload file {indir}/{id}-{kind}pld.parquet')
    sensor = pl.read_parquet(f'{indir}/{id}-{kind}pld.parquet')
    sensor = _remove_fill_values(sensor)

    # build a new data set based on info in `deploymentyaml.`
    # We will use ctd as the interpolant
    ds = xr.Dataset()
    attr = {}
    name = 'time'
    for atts in ncvar[name].keys():
        if atts != 'coordinates':
            attr[atts] = ncvar[name][atts]

    # If present, use the timebase specified in ncvar: timebase in the
    # deployment yaml.
    if 'timebase' not in ncvar:
        raise ValueError("Must specify timebase:source in netcdf_variables section of deployment yaml")
    if ncvar['timebase']['source'] not in sensor.columns:
        raise ValueError(f"timebase source: {ncvar['timebase']['source']} not found in pld1 columns")
    vals = sensor.select([ncvar['timebase']['source']]).to_numpy()[:, 0]
    indctd = np.where(~np.isnan(vals))[0]
    ds['time'] = (('time'), sensor.select('time').to_numpy()[indctd, 0].astype('datetime64[ns]'), attr)
    thenames = list(ncvar.keys())
    # Check yaml to see if interpolate has been set to True
    if "interpolate" in thenames:
        if ncvar["interpolate"]:
            interpolate = True
    for i in ['time', 'timebase', 'keep_variables', 'interpolate']:
        if i in thenames:
            thenames.remove(i)
    for name in thenames:
        _log.info('interpolating ' + name)
        if not ('method' in ncvar[name].keys()):
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
                    coarse_ints = np.arange(0, len(sensor) / coarsen_time, 1 / coarsen_time).astype(int)
                    sensor_sub = sensor.with_columns(pl.lit(coarse_ints).alias("coarse_ints"))
                    # Subsample the variable data keeping only the samples from the coarsened timeseries
                    sensor_sub_grouped = sensor_sub.with_columns(
                        pl.col('time').to_physical()
                    ).group_by(
                        pl.col('coarse_ints'),
                        maintain_order=True
                    ).mean().with_columns(
                        pl.col('time').cast(pl.Datetime('ms'))
                    )[:-1, :]
                    val2 = sensor_sub_grouped.select(sensorname).to_numpy()[:, 0]
                    val = _interp_gli_to_pld(sensor_sub_grouped, sensor, val2, indctd)
                if interpolate and not np.isnan(val).all():
                    time_original = sensor.select('time').to_numpy()[:, 0]
                    time_var = time_original[np.where(~np.isnan(val))[0]]
                    var_non_nan = val[np.where(~np.isnan(val))[0]]
                    time_timebase = sensor.select('time').to_numpy()[indctd, 0]
                    if val.dtype == "<M8[us]":
                        # for datetime, must convert to numerical, interpolate, then convert back
                        us_since_1970 = (var_non_nan - np.datetime64("1970-01-01")).astype(int)
                        val_int = np.interp(time_timebase.astype(float), time_var.astype(float), us_since_1970)
                        val_us = val_int.astype("timedelta64[us]")
                        val = np.datetime64("1970-01-01") + val_us
                    else:
                        val = np.interp(time_timebase.astype(float), time_var.astype(float), var_non_nan)

                    # interpolate only over those gaps that are smaller than 'maxgap'
                    # apparently maxgap is to be in somethng like seconds, and this data is in ms.  Certainly
                    # the default of 0.3 s was not intended.  Changing default to 10 s:
                    tg_ind = utils.find_gaps(time_var.astype(float), time_timebase.astype(float), maxgap*1000)
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

    # somehow this comes out unsorted:
    ds = ds.sortby(ds.time)
    # Drop duplicate timestamps and check how many are removed this way
    len_before_drop = len(ds.time)
    if hasattr(ds, "drop_duplicates"):
        ds = ds.drop_duplicates(dim="time")
    else:
        time_diffs = (ds.time.astype(int).diff(dim="time") > 1e-6).values
        time_diffs_list = list(time_diffs)
        time_diffs_list.append(True)
        good = np.array(time_diffs_list)
        ds = ds.isel(time=good)
    len_after_drop = len(ds.time)
    proportion_kept = len_after_drop / len_before_drop
    loss_str = f"{100 * (1 - proportion_kept)} % samples removed by timestamp deduplication."
    if proportion_kept < 0.5:
        raise ValueError(f"{loss_str} Check input data for duplicate timestamps")
    elif proportion_kept < 0.999:
        _log.warning(loss_str)
    else:
        _log.info(loss_str)

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
    ds.to_netcdf(outname, 'w',
                 encoding={'time': {'units': 'seconds since 1970-01-01T00:00:00Z',
                                    'dtype': 'float64'}})
    return outname


# aliases:
raw_to_L1timeseries = raw_to_L0timeseries = raw_to_timeseries
merge_rawnc = merge_parquet

__all__ = ['raw_to_rawnc', 'merge_parquet', 'raw_to_timeseries']
