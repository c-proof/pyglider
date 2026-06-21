"""
CDL/NC-based regression tests for the example-seaexplorer OG 1.0 timeseries.

Golden files live in tests/expected/ as cleaned .nc files with .cdl headers.
Tests compare pipeline output against golden NC using xarray.testing.assert_identical.

To regenerate golden files after an intentional change, run:
    python tests/_generate_expected_cdl.py
"""

from pathlib import Path

import numpy as np
import xarray as xr

import pyglider.seaexplorer as seaexplorer
import pyglider.ncprocess as ncprocess

from nc_test_helpers import assert_datasets_equal, compliance_msgs, run_compliance

library_dir = Path(__file__).parent.parent.absolute()
example_dir = library_dir / 'tests/example-data/'
expected_dir = library_dir / 'tests/expected/'

# ---------------------------------------------------------------------------
# Run the pipeline once at module level.
# raw_to_rawnc and merge_parquet are shared with the standard pipeline test;
# the parquet files in realtime_rawnc/ are reused here.
# ---------------------------------------------------------------------------
rawdir = str(example_dir / 'example-seaexplorer/realtime_raw/') + '/'
rawncdir = str(example_dir / 'example-seaexplorer/realtime_rawnc/') + '/'
deploymentyaml = str(example_dir / 'example-seaexplorer/deploymentRealtime_og10.yml')
l0tsdir = str(example_dir / 'example-seaexplorer/L0-timeseries-og10/') + '/'
griddir = str(example_dir / 'example-seaexplorer/L0-gridfiles-og10/') + '/'

seaexplorer.raw_to_rawnc(rawdir, rawncdir, deploymentyaml)
seaexplorer.merge_parquet(rawncdir, rawncdir, deploymentyaml, kind='sub')
outname = seaexplorer.raw_to_timeseries(rawncdir, l0tsdir, deploymentyaml, kind='sub')
outname_grid = ncprocess.make_gridfiles(outname, griddir, deploymentyaml)

output = xr.open_dataset(outname).load()
output_grid = xr.open_dataset(outname_grid).load()

expected_nc = xr.open_dataset(
    expected_dir / 'example-seaexplorer/L0-timeseries-og10/dfo-eva035-20190718.nc'
).load()
expected_grid = xr.open_dataset(
    expected_dir / 'example-seaexplorer/L0-gridfiles-og10/dfo-eva035-20190718.nc'
).load()


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_output_dimension():
    assert 'N_MEASUREMENTS' in output.dims, 'expected N_MEASUREMENTS dimension for OG 1.0'


def test_dataset():
    assert_datasets_equal(output, expected_nc)


def test_time_monotonic():
    """TIME must be strictly increasing."""
    t = output['TIME'].values.astype('datetime64[ns]').astype('float64')
    assert np.all(np.diff(t) > 0), 'TIME is not monotonically increasing'


def test_gridfile():
    assert_datasets_equal(output_grid, expected_grid)


def test_grid_time_monotonic():
    t = output_grid['time'].values.astype('datetime64[ns]').astype('float64')
    assert np.all(np.diff(t) > 0), 'grid time is not monotonically increasing'


# ---------------------------------------------------------------------------
# OG 1.0 compliance
# ---------------------------------------------------------------------------

def test_timeseries_og10_compliant():
    cc_data = run_compliance(outname, ['og'])
    result = cc_data['og']
    # it doesn't like that some pyglider variables have no `vocabulary` attribute; ignore
    # Variables should have valid vocabulary URIs: variable PROFILE_NUMBER should have attribute 'vocabulary', value is a URI from the OG1 collection
    assert result['high_count'] == 1, \
        "og high priority errors:\n" + "\n".join(compliance_msgs(result, 3))
    assert result['medium_count'] == 0, \
        "og medium priority errors:\n" + "\n".join(compliance_msgs(result, 2))

    cc_data = run_compliance(outname, ['cf:1.8'])
    result = cc_data['cf:1.8']
    assert result['high_count'] == 0, \
        "cf:1.8 high priority errors:\n" + "\n".join(compliance_msgs(result, 3))
    assert result['medium_count'] == 0, \
        "cf:1.8 medium priority errors:\n" + "\n".join(compliance_msgs(result, 2))
