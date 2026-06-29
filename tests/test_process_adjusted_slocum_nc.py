"""
CDL/NC-based regression tests for the slocum process_adjusted pipeline.

Golden files live in tests/expected/ as cleaned .nc files with .cdl headers.
Pipeline output is written to example-data/ (not expected/) so the golden
files are not overwritten during the test run.

To regenerate golden files after an intentional change, run:
    python tests/_generate_expected_cdl.py
"""

from pathlib import Path

import numpy as np
import xarray as xr

import pyglider.slocum as slocum
from pyglider.process_adjusted import run_process_adjusted

from nc_test_helpers import assert_datasets_equal

library_dir = Path(__file__).parent.parent.absolute()
expected_dir = library_dir / 'tests/expected/example-slocum'
yaml_dir = library_dir / 'tests/example-data/example-slocum'

# ---------------------------------------------------------------------------
# Run the pipeline once at module level.
# Write output to example-data/ so the golden NC in expected/ is not touched.
# ---------------------------------------------------------------------------
slocum.binary_to_timeseries(
    str(yaml_dir / 'realtime_raw/') + '/',
    yaml_dir / 'cac/',
    str(yaml_dir / 'L0-timeseries/') + '/',
    str(yaml_dir / 'deploymentRealtime.yml'),
    search='*.[s|t]bd',
    profile_filt_time=20,
    profile_min_time=20,
)

_output_base = yaml_dir / 'adjusted_output'
outname = run_process_adjusted(
    _output_base,
    deploy_name='dfo-rosie713-20190615',
    deployfile=yaml_dir / 'deploymentRealtime.yml',
    adjustedyaml=yaml_dir / 'adjusted.yml',
    input_dir=yaml_dir,
)
output = xr.open_dataset(outname).load()

gridfile = _output_base / 'L0-gridfiles' / 'dfo-rosie713-20190615_grid_adjusted.nc'
output_grid = xr.open_dataset(gridfile).load()

# ---------------------------------------------------------------------------
# Load golden NC files
# ---------------------------------------------------------------------------
expected = xr.open_dataset(
    expected_dir / 'L0-timeseries' / 'dfo-rosie713-20190615_adjusted.nc'
).load()
expected_grid = xr.open_dataset(
    expected_dir / 'L0-gridfiles' / 'dfo-rosie713-20190615_grid_adjusted.nc'
).load()


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_dataset():
    assert_datasets_equal(output, expected)


def test_dataset_grid():
    assert_datasets_equal(output_grid, expected_grid)


def test_grid_time_monotonic():
    t = output_grid['time'].values.astype('datetime64[ns]').astype('float64')
    assert np.all(np.diff(t) > 0), 'grid time is not monotonically increasing'
