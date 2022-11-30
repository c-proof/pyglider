import polars as pl
import pytest
from pathlib import Path
import os
library_dir = Path(__file__).parent.parent.absolute()
example_dir = str(library_dir) + '/tests/example-data/'
data_dir = str(library_dir) + '/tests/data/'
os.system('rm ' + data_dir + '/realtime_rawnc/*')

import pyglider.seaexplorer as seaexplorer

def test__outputname():
    fnout, filenum = seaexplorer._outputname(data_dir + '/realtime_raw/sea035.12.pld1.sub.36',
                                             data_dir + '/realtime_rawnc/')
    assert fnout == data_dir + '/realtime_rawnc/sea035.0012.pld1.sub.0036.parquet'
    assert filenum == 36


def test_raw_to_rawnc():
    # Default settings on a clean folder
    result_default = seaexplorer.raw_to_rawnc(data_dir + '/realtime_raw/',
                                              data_dir + '/realtime_rawnc/', None)
    # Test the reprocess flag works
    result_reprocess = seaexplorer.raw_to_rawnc(data_dir + '/realtime_raw/',
                                                data_dir + '/realtime_rawnc/',
                                                None, incremental=False)
    # Check that reprocessing not preformed by default
    result_no_new_files = seaexplorer.raw_to_rawnc(data_dir + '/realtime_raw/',
                                                   data_dir + '/realtime_rawnc/',
                                                   None)
    # Reject all payload files with fewer than 10000 lines
    result_strict = seaexplorer.raw_to_rawnc(data_dir + '/realtime_raw/',
                                             data_dir + '/realtime_rawnc/',
                                             None,
                                             incremental=False,
                                             min_samples_in_file=10000)
    assert result_default is True
    assert result_reprocess is True
    assert result_no_new_files is False
    assert result_strict is False


def test__needsupdating():
    ftype = 'pld1'
    fin = data_dir + '/realtime_raw/sea035.12.pld1.sub.36'
    fout = data_dir+ '/realtime_rawnc/sea035.0012.pld1.sub.0036.parquet'
    result_badpath = seaexplorer._needsupdating(ftype, fin, 'baz')
    result_goodpath = seaexplorer._needsupdating(ftype, fin, fout)
    assert result_badpath is True
    assert result_goodpath is False


def test_merge_rawnc():
    result_default = seaexplorer.merge_parquet(
            data_dir + '/realtime_rawnc/',
            data_dir + '/realtime_rawnc/',
            example_dir +  '/example-seaexplorer/deploymentRealtime.yml')

    result_sub = seaexplorer.merge_parquet(
            data_dir + '/realtime_rawnc/',
            data_dir + '/realtime_rawnc/',
            example_dir + '/example-seaexplorer/deploymentRealtime.yml',
            kind='sub')
    assert result_default is False
    assert result_sub is True


def test__interp_gli_to_pld():
    # function should interpolate values from the glider dataset to sampling frequency of payload dataset
    glider = pl.read_parquet(data_dir + '/realtime_rawnc/sea035.0012.gli.sub.0036.parquet')
    ds = pl.read_parquet(data_dir + '/realtime_rawnc/sea035.0012.pld1.sub.0036.parquet')
    val = glider.select("Pitch").to_numpy()[:, 0]
    pitch_interp = seaexplorer._interp_gli_to_pld(glider, ds, val, None)
    assert len(pitch_interp) == ds.shape[0]


def test_raw_to_timeseries():
    # Test default, will fail as we have sub data, not raw data
    with pytest.raises(FileNotFoundError) as missing_file_exc:
        result_default = seaexplorer.raw_to_timeseries(data_dir + '/realtime_rawnc/',
                                                        data_dir + '/l0-profiles/',
                                                        example_dir + '/example-seaexplorer/deploymentRealtime.yml',
                                                        )
    result_sub = seaexplorer.raw_to_timeseries(data_dir + '/realtime_rawnc/',
                                                    data_dir + '/l0-profiles/',
                                                    example_dir + '/example-seaexplorer/deploymentRealtime.yml',
                                                    kind='sub')
    assert 'No such file or directory' in str(missing_file_exc)
    assert result_sub == data_dir + '/l0-profiles/dfo-eva035-20190718.nc'

