"""
CDL/NC-based regression tests for the example-slocum L0 timeseries.

Golden files live in tests/expected/ as cleaned .nc files with .cdl headers.
Tests compare pipeline output against golden NC using xarray.testing.assert_identical.

To regenerate golden files after an intentional change, run:
    python tests/_generate_expected_cdl.py
"""

from pathlib import Path

import numpy as np
import xarray as xr

import pyglider.ncprocess as ncprocess
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
deploymentyaml = str(example_dir / 'example-slocum/deploymentRealtime.yml')
tsdir = str(example_dir / 'example-slocum/L0-timeseries/') + '/'
profiledir = str(example_dir / 'example-slocum/L0-profiles/')

outname = slocum.binary_to_timeseries(
    binarydir,
    cacdir,
    tsdir,
    deploymentyaml,
    search='*.[s|t]bd',
    profile_filt_time=20,
    profile_min_time=20,
)
ncprocess.extract_timeseries_profiles(outname, profiledir, deploymentyaml, force=True)

output = xr.open_dataset(outname).load()

expected_nc = xr.open_dataset(
    expected_dir / 'example-slocum/L0-timeseries/dfo-rosie713-20190615.nc'
).load()


# ---------------------------------------------------------------------------
# Timeseries tests
# ---------------------------------------------------------------------------

def test_dataset():
    assert_datasets_equal(output, expected_nc)


def test_time_monotonic():
    """Time must be strictly increasing."""
    t = output['time'].values.astype('datetime64[ns]').astype('float64')
    assert np.all(np.diff(t) > 0), 'time is not monotonically increasing'


# ---------------------------------------------------------------------------
# Compliance tests (ported from test_slocum.py)
# ---------------------------------------------------------------------------

def test_profiles_compliant():
    path = Path(profiledir) / 'dfo-rosie713-20190620T1313.nc'
    cc_data = run_compliance(path, ['gliderdac', 'cf:1.8'])
    for checker in ['gliderdac', 'cf:1.8']:
        result = cc_data[checker]
        assert result['high_count'] == 0, \
            f"{checker} high priority errors:\n" + "\n".join(compliance_msgs(result, 3))
        assert result['medium_count'] == 0, \
            f"{checker} medium priority errors:\n" + "\n".join(compliance_msgs(result, 2))
        assert result['low_count'] == 0, \
            f"{checker} low priority errors:\n" + "\n".join(compliance_msgs(result, 1))


def test_timeseries_compliant():
    cc_data = run_compliance(outname, ['cf:1.8'])
    result = cc_data['cf:1.8']
    assert result['high_count'] == 0, \
        "cf:1.8 high priority errors:\n" + "\n".join(compliance_msgs(result, 3))
    assert result['medium_count'] == 0, \
        "cf:1.8 medium priority errors:\n" + "\n".join(compliance_msgs(result, 2))
    assert result['low_count'] == 0, \
        "cf:1.8 low priority errors:\n" + "\n".join(compliance_msgs(result, 1))
