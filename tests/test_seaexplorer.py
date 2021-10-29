import xarray as xr
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
