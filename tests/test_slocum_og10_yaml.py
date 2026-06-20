"""
YAML-summary-based regression tests for the example-slocum OG 1.0 timeseries.

Tests cover:
  - variable list
  - global attributes (excluding dynamic/version fields)
  - per-variable CF attributes
  - per-variable statistics (min, max, mean, std, n_valid, n_nan)
  - per-variable time-derivative statistics (diff_mean, diff_std)
  - TIME monotonicity
  - OG 1.0 compliance (cc-plugin-og)

To regenerate the YAML after an intentional change, run:
    python tests/_generate_expected_yaml.py
"""

from pathlib import Path
from unittest import result

import numpy as np
import pytest
import xarray as xr
import yaml

import pyglider.slocum as slocum

from yaml_test_helpers import ATOL, RTOL, SKIP_ATTRS, compliance_msgs, run_compliance, stats, to_float

library_dir = Path(__file__).parent.parent.absolute()
example_dir = library_dir / 'tests/example-data/'

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

output = xr.open_dataset(outname)

# ---------------------------------------------------------------------------
# Load expected YAML summary
# ---------------------------------------------------------------------------
expected_yaml_path = (
    library_dir
    / 'tests/expected/example-slocum/L0-timeseries-og10/dfo-rosie713-20190615.yml'
)
with open(expected_yaml_path) as f:
    expected = yaml.safe_load(f)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_output_dimension():
    assert 'N_MEASUREMENTS' in output.dims, 'expected N_MEASUREMENTS dimension for OG 1.0'


def test_variable_list():
    assert sorted(output.variables) == sorted(expected['variables'])


def test_global_attrs():
    out_attrs = {k: str(v) for k, v in output.attrs.items() if k not in SKIP_ATTRS}
    exp_attrs = {k: v for k, v in expected['attrs'].items() if k not in SKIP_ATTRS}
    assert out_attrs == exp_attrs


@pytest.mark.parametrize('var', sorted(expected['variables']))
def test_variable_attrs(var):
    exp_attrs = expected['variables'][var]['attrs']
    out_attrs = {k: str(v) for k, v in output[var].attrs.items()}
    assert out_attrs == exp_attrs


@pytest.mark.parametrize('var', sorted(expected['variables']))
def test_variable_stats(var):
    vals = to_float(output[var])
    actual = stats(vals)
    exp = expected['variables'][var]

    assert actual['n_valid'] == exp['n_valid'], f'{var}: n_valid mismatch'
    assert actual['n_nan'] == exp['n_nan'], f'{var}: n_nan mismatch'

    for stat in ('min', 'max', 'mean', 'std'):
        if stat not in exp:
            continue
        np.testing.assert_allclose(
            actual[stat], exp[stat], rtol=RTOL, atol=ATOL,
            err_msg=f'{var}: {stat} mismatch',
        )


@pytest.mark.parametrize('var', sorted(expected['variables']))
def test_variable_diffstats(var):
    """Time-derivative statistics catch temporal scrambling."""
    vals = to_float(output[var])
    actual = stats(vals)
    exp = expected['variables'][var]

    for stat in ('diff_mean', 'diff_std'):
        if stat not in exp:
            continue
        np.testing.assert_allclose(
            actual[stat], exp[stat], rtol=RTOL, atol=ATOL,
            err_msg=f'{var}: {stat} mismatch',
        )


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
