"""
CDL/NC-based regression tests for the seaexplorer process_adjusted pipeline.

Golden files live in tests/expected/ as cleaned .nc files with .cdl headers.
Pipeline output is written to example-data/ (not expected/) so the golden
files are not overwritten during the test run.

To regenerate golden files after an intentional change, run:
    python tests/_generate_expected_cdl.py
"""

from pathlib import Path

import xarray as xr

import pyglider.seaexplorer as seaexplorer
from pyglider.process_adjusted import run_process_adjusted

from nc_test_helpers import assert_datasets_equal

library_dir = Path(__file__).parent.parent.absolute()
expected_dir = library_dir / 'tests/expected/example-seaexplorer'
yaml_dir = library_dir / 'tests/example-data/example-seaexplorer'

# ---------------------------------------------------------------------------
# Run the pipeline once at module level.
# Write output to example-data/ so the golden NC in expected/ is not touched.
# ---------------------------------------------------------------------------
_rawdir = str(yaml_dir / 'realtime_raw/') + '/'
_rawncdir = str(yaml_dir / 'realtime_rawnc/') + '/'
_deploymentyaml = str(yaml_dir / 'deploymentRealtime.yml')
_l0tsdir = str(yaml_dir / 'L0-timeseries/') + '/'
seaexplorer.raw_to_rawnc(_rawdir, _rawncdir, _deploymentyaml)
seaexplorer.merge_parquet(_rawncdir, _rawncdir, _deploymentyaml, kind='sub')
seaexplorer.raw_to_timeseries(_rawncdir, _l0tsdir, _deploymentyaml, kind='sub')

# Write adjusted outputs to a subdirectory of example-data, not to expected/
_output_base = yaml_dir / 'adjusted_output'
outname = run_process_adjusted(
    _output_base,
    deploy_name='dfo-eva035-20190718',
    deployfile=yaml_dir / 'deploymentRealtime.yml',
    adjustedyaml=yaml_dir / 'adjusted.yml',
    input_dir=yaml_dir,
)
output = xr.open_dataset(outname).load()

gridfile = _output_base / 'L0-gridfiles' / 'dfo-eva035-20190718_grid_adjusted.nc'
output_grid = xr.open_dataset(gridfile).load()

# ---------------------------------------------------------------------------
# Load golden NC files
# ---------------------------------------------------------------------------
expected = xr.open_dataset(
    expected_dir / 'L0-timeseries' / 'dfo-eva035-20190718_adjusted.nc'
).load()
expected_grid = xr.open_dataset(
    expected_dir / 'L0-gridfiles' / 'dfo-eva035-20190718_grid_adjusted.nc'
).load()


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_dataset():
    assert_datasets_equal(output, expected)


def test_dataset_grid():
    assert_datasets_equal(output_grid, expected_grid)


def test_grid_time_monotonic():
    import numpy as np
    t = output_grid['time'].values.astype('datetime64[ns]').astype('float64')
    assert np.all(np.diff(t) > 0), 'grid time is not monotonically increasing'
