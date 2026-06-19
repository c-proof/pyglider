"""
YAML-summary-based regression tests for the slocum process_adjusted pipeline.

Intended to replace test_process_adjusted_slocum.py.  Tests cover:
  - variable list
  - global attributes (excluding dynamic fields)
  - per-variable CF attributes
  - per-variable statistics (min, max, mean, std, n_valid, n_nan)
  - per-variable time-derivative statistics (diff_mean, diff_std)
  - L0-gridfiles output from the same pipeline

To regenerate the YAML golden files after an intentional change, run:
    python tests/_generate_expected_yaml.py
"""

from pathlib import Path

import numpy as np
import pytest
import xarray as xr
import yaml

from pyglider.process_adjusted import run_process_adjusted

library_dir = Path(__file__).parent.parent.absolute()
expected_dir = library_dir / 'tests/expected/example-slocum'
yaml_dir = library_dir / 'tests/example-data/example-slocum'

# ---------------------------------------------------------------------------
# Run the pipeline once at module level
# ---------------------------------------------------------------------------
outname = run_process_adjusted(
    expected_dir,
    deploy_name='dfo-rosie713-20190615',
    deployfile=yaml_dir / 'deploymentRealtime.yml',
    adjustedyaml=yaml_dir / 'adjusted.yml',
    input_dir=yaml_dir,
)
output = xr.open_dataset(outname)

# L0-gridfile produced by the same pipeline call
gridfile = expected_dir / 'L0-gridfiles' / 'dfo-rosie713-20190615_grid_adjusted.nc'
output_grid = xr.open_dataset(gridfile)

# ---------------------------------------------------------------------------
# Load expected YAML summaries
# ---------------------------------------------------------------------------
expected_yaml_path = (
    expected_dir / 'L0-timeseries' / 'dfo-rosie713-20190615_adjusted.yml'
)
with open(expected_yaml_path) as f:
    expected = yaml.safe_load(f)

expected_grid_yaml_path = (
    expected_dir / 'L0-gridfiles' / 'dfo-rosie713-20190615_grid_adjusted.yml'
)
with open(expected_grid_yaml_path) as f:
    expected_grid = yaml.safe_load(f)

from yaml_test_helpers import RTOL, ATOL, SKIP_ATTRS, to_float, stats


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_variable_list():
    assert sorted(output.variables) == sorted(expected['variables'])


def test_global_attrs():
    out_attrs = {k: str(v) for k, v in output.attrs.items() if k not in SKIP_ATTRS}
    exp_attrs = {k: v for k, v in expected['attrs'].items() if k not in SKIP_ATTRS}
    assert out_attrs == exp_attrs


@pytest.mark.parametrize('var', sorted(expected['variables']))
def test_variable_attrs(var):
    assert {k: str(v) for k, v in output[var].attrs.items()} == expected['variables'][var]['attrs']


@pytest.mark.parametrize('var', sorted(expected['variables']))
def test_variablestats(var):
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


# ---------------------------------------------------------------------------
# L0-gridfile tests
# ---------------------------------------------------------------------------

def test_grid_variable_list():
    assert sorted(output_grid.variables) == sorted(expected_grid['variables'])


def test_grid_global_attrs():
    out_attrs = {k: str(v) for k, v in output_grid.attrs.items() if k not in SKIP_ATTRS}
    exp_attrs = {k: v for k, v in expected_grid['attrs'].items() if k not in SKIP_ATTRS}
    assert out_attrs == exp_attrs


@pytest.mark.parametrize('var', sorted(expected_grid['variables']))
def test_grid_variable_attrs(var):
    assert {k: str(v) for k, v in output_grid[var].attrs.items()} == expected_grid['variables'][var]['attrs']


@pytest.mark.parametrize('var', sorted(expected_grid['variables']))
def test_grid_variablestats(var):
    vals = to_float(output_grid[var])
    actual = stats(vals)
    exp = expected_grid['variables'][var]
    assert actual['n_valid'] == exp['n_valid'], f'{var}: n_valid mismatch'
    assert actual['n_nan'] == exp['n_nan'], f'{var}: n_nan mismatch'
    for stat in ('min', 'max', 'mean', 'std'):
        if stat not in exp:
            continue
        np.testing.assert_allclose(
            actual[stat], exp[stat], rtol=RTOL, atol=ATOL,
            err_msg=f'{var}: {stat} mismatch',
        )


@pytest.mark.parametrize('var', sorted(expected_grid['variables']))
def test_grid_variable_diffstats(var):
    vals = to_float(output_grid[var])
    actual = stats(vals)
    exp = expected_grid['variables'][var]
    for stat in ('diff_mean', 'diff_std'):
        if stat not in exp:
            continue
        np.testing.assert_allclose(
            actual[stat], exp[stat], rtol=RTOL, atol=ATOL,
            err_msg=f'{var}: {stat} mismatch',
        )


def test_grid_time_monotonic():
    t = output_grid['time'].values.astype('datetime64[ns]').astype('float64')
    assert np.all(np.diff(t) > 0), 'grid time is not monotonically increasing'
