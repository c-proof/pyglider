"""
CDL/NC-based regression tests for the example-slocum OG 1.0 timeseries.

Golden files live in tests/expected/ as cleaned .nc files with .cdl headers.
Tests compare pipeline output against golden NC using xarray.testing.assert_identical.

To regenerate golden files after an intentional change, run:
    python tests/_generate_expected_cdl.py
"""

from pathlib import Path

import numpy as np
import xarray as xr

import pyglider.slocum as slocum

from nc_test_helpers import assert_datasets_equal, compliance_msgs, run_compliance

library_dir = Path(__file__).parent.parent.absolute()
example_dir = library_dir / 'tests/example-data/'
expected_dir = library_dir / 'tests/expected/'

# ---------------------------------------------------------------------------
# Run the pipeline once at module level
# ---------------------------------------------------------------------------
cacdir = example_dir / 'example-slocum/cac/'
binarydir = str(example_dir / 'example-slocum/realtime_raw/') + '/'
deploymentyaml = str(example_dir / 'example-slocum/deploymentRealtime_og10.yml')
tsdir = str(example_dir / 'example-slocum/L0-timeseries-og10/') + '/'

outname = slocum.binary_to_timeseries(
    binarydir,
    cacdir,
    tsdir,
    deploymentyaml,
    search='*.[s|t]bd',
    profile_filt_time=20,
    profile_min_time=20,
)

output = xr.open_dataset(outname).load()

expected_nc = xr.open_dataset(
    expected_dir / 'example-slocum/L0-timeseries-og10/dfo-rosie713-20190615.nc'
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


# ---------------------------------------------------------------------------
# OG 1.0 compliance
# ---------------------------------------------------------------------------

def test_timeseries_og10_compliant():
    cc_data = run_compliance(outname, ['og'])
    result = cc_data['og']
    # it doesn't like that some pyglider variables have no `vocabulary` attribute; ignore
    # Variables should have valid vocabulary URIs: variable WAYPOINT_LATITUDE should have attribute 'vocabulary', value is a URI from the OG1 collection
    assert result['high_count'] == 1, \
        "og high priority errors:\n" + "\n".join(compliance_msgs(result, 3))
    assert result['medium_count'] == 0, \
        "og medium priority errors:\n" + "\n".join(compliance_msgs(result, 2))
    assert result['low_count'] == 0, \
        "og low priority errors:\n" + "\n".join(compliance_msgs(result, 1))

    cc_data = run_compliance(outname, ['cf:1.8'])
    result = cc_data['cf:1.8']
    assert result['high_count'] == 0, \
        "cf:1.8 high priority errors:\n" + "\n".join(compliance_msgs(result, 3))
    assert result['medium_count'] == 0, \
        "cf:1.8 medium priority errors:\n" + "\n".join(compliance_msgs(result, 2))
