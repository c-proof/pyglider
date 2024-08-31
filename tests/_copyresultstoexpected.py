"""
If we run the tests and decide we want all the test results to be copied over, this is faster...
"""
from pathlib import Path
from shutil import copy

todo = {'example-data/example-seaexplorer/L0-timeseries-test/dfo-eva035-20190718.nc':
        'expected/example-seaexplorer/L0-timeseries',
        'example-data/example-seaexplorer-raw/L0-timeseries-test/dfo-bb046-20200908.nc': 'expected/example-seaexplorer-raw/L0-timeseries',
        'example-data/example-slocum/L0-timeseries/dfo-rosie713-20190615.nc': 'expected/example-slocum/L0-timeseries',
        'example-data/example-slocum-littleendian/L0-timeseries-test/dfo-maria997-20220614.nc': 'expected/example-slocum-littleendian/L0-timeseries'

        }

for td in todo:
    copy(td, todo[td])

