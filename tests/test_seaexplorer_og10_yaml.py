"""
YAML-summary-based regression tests for the example-seaexplorer OG 1.0 timeseries.

Tests cover:
  - variable list
  - global attributes (excluding dynamic/version fields)
  - per-variable CF attributes
  - per-variable statistics (min, max, mean, std, n_valid, n_nan)
  - per-variable time-derivative statistics (diff_mean, diff_std)
  - TIME monotonicity

To regenerate the YAML after an intentional change, run:
    python tests/_generate_expected_yaml.py
"""

from pathlib import Path

import numpy as np
import pytest
import xarray as xr
import yaml

import pyglider.seaexplorer as seaexplorer

from yaml_test_helpers import ATOL, RTOL, SKIP_ATTRS, stats, to_float

library_dir = Path(__file__).parent.parent.absolute()
example_dir = library_dir / 'tests/example-data/'

# ---------------------------------------------------------------------------
# Run the pipeline once at module level.
# raw_to_rawnc and merge_parquet are shared with the standard pipeline test;
# the parquet files in realtime_rawnc/ are reused here.
# ---------------------------------------------------------------------------
rawdir = str(example_dir / 'example-seaexplorer/realtime_raw/') + '/'
rawncdir = str(example_dir / 'example-seaexplorer/realtime_rawnc/') + '/'
deploymentyaml = str(example_dir / 'example-seaexplorer/deploymentRealtime_og10.yml')
l0tsdir = str(example_dir / 'example-seaexplorer/L0-timeseries-og10/') + '/'

seaexplorer.raw_to_rawnc(rawdir, rawncdir, deploymentyaml)
seaexplorer.merge_parquet(rawncdir, rawncdir, deploymentyaml, kind='sub')
outname = seaexplorer.raw_to_timeseries(rawncdir, l0tsdir, deploymentyaml, kind='sub')

output = xr.open_dataset(outname)

# ---------------------------------------------------------------------------
# Load expected YAML summary
# ---------------------------------------------------------------------------
expected_yaml_path = (
    library_dir
    / 'tests/expected/example-seaexplorer/L0-timeseries-og10/dfo-eva035-20190718.yml'
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
