import xarray as xr
import logging
from pathlib import Path
import sys

library_dir = Path(__file__).parent.parent.absolute()
sys.path.append(str(library_dir))
import pyglider.seaexplorer as seaexplorer


def test__outputname():
    fnout, filenum = seaexplorer._outputname('example-seaexplorer/realtime_raw/sea035.12.pld1.sub.36',
                            'example-seaexplorer/realtime_rawnc/')
    assert fnout == 'example-seaexplorer/realtime_rawnc/sea035.0012.pld1.sub.0036.nc'
    assert filenum == 36
   
   
def test__needsupdating():
    ftype = 'pld1'
    fin = 'example-seaexplorer/realtime_raw/sea035.12.pld1.sub.36'
    fout = 'example-seaexplorer/realtime_rawnc/sea035.0012.pld1.sub.0036.nc'
    result_badpath = seaexplorer._needsupdating(ftype, fin, 'baz')
    result_goodpath = seaexplorer._needsupdating(ftype, fin, fout)
    assert result_badpath is True
    assert result_goodpath is False


def test_raw_to_rawnc():
    # Default settings on a clean folder
    result_default = seaexplorer.raw_to_rawnc('example-seaexplorer/realtime_raw/',
                             'example-seaexplorer/realtime_rawnc/',
                             None)
    # Test the reprocess flag works
    result_reprocess = seaexplorer.raw_to_rawnc('example-seaexplorer/realtime_raw/',
                                              'example-seaexplorer/realtime_rawnc/',
                                              None,
                                              incremental=False)
    # Check that reprocessing not preformed by default
    result_no_new_files = seaexplorer.raw_to_rawnc('example-seaexplorer/realtime_raw/',
                             'example-seaexplorer/realtime_rawnc/',
                             None)
    # Reject all payload files with fewer than 10000 lines
    result_strict = seaexplorer.raw_to_rawnc('example-seaexplorer/realtime_raw/',
                                              'example-seaexplorer/realtime_rawnc/',
                                             None,
                                             incremental=False,
                                             min_samples_in_file=10000)
    assert result_default is True
    assert result_reprocess is True
    assert result_no_new_files is False
    assert result_strict is False


def test_merge_rawnc():
    result_default = seaexplorer.merge_rawnc('example-seaexplorer/realtime_rawnc/',
                                             'example-seaexplorer/realtime_rawnc/',
                                             'example-seaexplorer/deploymentRealtime.yml')

    result_sub = seaexplorer.merge_rawnc('example-seaexplorer/realtime_rawnc/',
                                             'example-seaexplorer/realtime_rawnc/',
                                             'example-seaexplorer/deploymentRealtime.yml',
                                             kind='sub')
    assert result_default is False
    assert result_sub is True


def test__interp_gli_to_pld():
    pass


def test__interp_pld_to_pld():
    pass


def test_raw_to_L0timeseries():
    pass
    
    