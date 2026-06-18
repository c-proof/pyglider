import os
from pathlib import Path

import numpy as np
import polars as pl
import pytest
import yaml

import pyglider.seaexplorer as seaexplorer

os.system('rm tests/data/realtime_rawnc/*')
library_dir = Path(__file__).parent.parent.absolute()
example_dir = library_dir / 'tests/example-data/'


def test__outputname():
    fnout, filenum = seaexplorer._outputname(
        'tests/data/realtime_raw/sea035.12.pld1.sub.36', 'tests/data/realtime_rawnc/'
    )
    assert fnout == 'tests/data/realtime_rawnc/sea035.0012.pld1.sub.0036.parquet'
    assert filenum == 36


def test_raw_to_rawnc():
    # Default settings on a clean folder
    result_default = seaexplorer.raw_to_rawnc(
        'tests/data/realtime_raw/', 'tests/data/realtime_rawnc/', None
    )
    assert result_default is True
    # Test the reprocess flag works
    result_reprocess = seaexplorer.raw_to_rawnc(
        'tests/data/realtime_raw/',
        'tests/data/realtime_rawnc/',
        None,
        incremental=False,
    )
    assert result_reprocess is True

    # Check that reprocessing not preformed by default
    result_no_new_files = seaexplorer.raw_to_rawnc(
        'tests/data/realtime_raw/', 'tests/data/realtime_rawnc/', None
    )
    assert result_no_new_files is False

    # Reject all payload files with fewer than 10000 lines
    result_strict = seaexplorer.raw_to_rawnc(
        'tests/data/realtime_raw/',
        'tests/data/realtime_rawnc/',
        None,
        incremental=False,
        min_samples_in_file=10000,
    )
    assert result_strict is False


def test__needsupdating():
    ftype = 'pld1'
    fin = 'tests/data/realtime_raw/sea035.12.pld1.sub.36'
    fout = 'tests/data/realtime_rawnc/sea035.0012.pld1.sub.0036.parquet'
    result_badpath = seaexplorer._needsupdating(ftype, fin, 'baz')
    result_goodpath = seaexplorer._needsupdating(ftype, fin, fout)
    assert result_badpath is True
    assert result_goodpath is False


def test_merge_rawnc():
    result_default = seaexplorer.merge_parquet(
        'tests/data/realtime_rawnc/',
        'tests/data/realtime_rawnc/',
        str(example_dir / 'example-seaexplorer/deploymentRealtime.yml'),
    )

    result_sub = seaexplorer.merge_parquet(
        'tests/data/realtime_rawnc/',
        'tests/data/realtime_rawnc/',
        str(example_dir / 'example-seaexplorer/deploymentRealtime.yml'),
        kind='sub',
    )
    assert result_default is False
    assert result_sub is True


def test__remove_fill_values():
    # This should convert values equallling 9999 in the original df to nan
    df_in = pl.read_parquet(
        'tests/data/realtime_rawnc/sea035.0012.pld1.sub.0036.parquet'
    )
    df_out = seaexplorer._remove_fill_values(df_in)
    assert (df_in.select('GPCTD_DOF').to_numpy()[:, 0] == 9999).all()
    assert np.isnan(df_out.select('GPCTD_DOF').to_numpy()[:, 0]).all()


def test__interp_gli_to_pld():
    # function should interpolate values from the glider dataset to sampling frequency of payload dataset
    glider = pl.read_parquet(
        'tests/data/realtime_rawnc/sea035.0012.gli.sub.0036.parquet'
    )
    ds = pl.read_parquet('tests/data/realtime_rawnc/sea035.0012.pld1.sub.0036.parquet')
    val = glider.select('Pitch').to_numpy()[:, 0]
    pitch_interp = seaexplorer._interp_gli_to_pld(glider, ds, val, None)
    assert len(pitch_interp) == ds.shape[0]


def test_raw_to_timeseries():
    # Test default, will fail as we have sub data, not raw data
    with pytest.raises(FileNotFoundError) as missing_file_exc:
        result_default = seaexplorer.raw_to_timeseries(
            'tests/data/realtime_rawnc/',
            'tests/data/l0-profiles/',
            str(example_dir / 'example-seaexplorer/deploymentRealtime.yml'),
        )
    result_sub = seaexplorer.raw_to_timeseries(
        'tests/data/realtime_rawnc/',
        'tests/data/l0-profiles/',
        str(example_dir / 'example-seaexplorer/deploymentRealtime.yml'),
        kind='sub',
    )
    assert 'No such file or directory' in str(missing_file_exc)
    assert result_sub == 'tests/data/l0-profiles/dfo-eva035-20190718.nc'


def test_missing_bad_timebase():
    # Prepare yaml files with bad timebase and no timebase
    with open(example_dir / 'example-seaexplorer/deploymentRealtime.yml') as fin:
        deployment = yaml.safe_load(fin)
    deployment['netcdf_variables']['timebase']['source'] = 'non existing sensor'
    with open(example_dir / 'example-seaexplorer/bad_timebase.yml', 'w') as fin:
        yaml.dump(deployment, fin)
    deployment['netcdf_variables'].pop('timebase')
    with open(example_dir / 'example-seaexplorer/no_timebase.yml', 'w') as fin:
        yaml.dump(deployment, fin)
    with pytest.raises(ValueError) as bad_timebase_exc:
        result_bad_timebase = seaexplorer.raw_to_timeseries(
            'tests/data/realtime_rawnc/',
            'tests/data/l0-profiles/',
            str(example_dir / 'example-seaexplorer/bad_timebase.yml'),
            kind='sub',
        )
    with pytest.raises(ValueError) as no_timebase_exc:
        result_no_timebase = seaexplorer.raw_to_timeseries(
            'tests/data/realtime_rawnc/',
            'tests/data/l0-profiles/',
            str(example_dir / 'example-seaexplorer/no_timebase.yml'),
            kind='sub',
        )
    assert 'sensor not found in pld1 columns' in str(bad_timebase_exc)
    assert 'Must specify timebase' in str(no_timebase_exc)
