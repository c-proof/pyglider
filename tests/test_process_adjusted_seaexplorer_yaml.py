"""
YAML-summary-based regression tests for the seaexplorer process_adjusted pipeline.

Intended to replace test_process_adjusted_seaexplorer.py.  Tests cover:
  - variable list
  - global attributes (excluding dynamic fields)
  - per-variable CF attributes
  - per-variable statistics (min, max, mean, std, n_valid, n_nan)
  - per-variable time-derivative statistics (diff_mean, diff_std)

To regenerate the YAML golden file after an intentional change, run:
    python tests/_generate_expected_yaml.py
"""

from pathlib import Path

import numpy as np
import pytest
import xarray as xr
import yaml

from pyglider.process_adjusted import run_process_adjusted

library_dir = Path(__file__).parent.parent.absolute()
expected_dir = library_dir / 'tests/expected/example-seaexplorer'
yaml_dir = library_dir / 'tests/example-data/example-seaexplorer'

# ---------------------------------------------------------------------------
# Run the pipeline once at module level
# ---------------------------------------------------------------------------
outname = run_process_adjusted(
    expected_dir,
    deploy_name='dfo-eva035-20190718',
    deployfile=yaml_dir / 'deploymentRealtime.yml',
    adjustedyaml=yaml_dir / 'adjusted.yml',
)
output = xr.open_dataset(outname)

# ---------------------------------------------------------------------------
# Load expected YAML summary
# ---------------------------------------------------------------------------
expected_yaml_path = (
    expected_dir / 'L0-timeseries' / 'dfo-eva035-20190718_adjusted.yml'
)
with open(expected_yaml_path) as f:
    expected = yaml.safe_load(f)

RTOL = 1e-5
ATOL = 0.0
SKIP_ATTRS = {'date_created', 'date_issued', 'history'}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _to_float(da):
    vals = da.values
    if np.issubdtype(vals.dtype, np.datetime64):
        vals = vals.astype('datetime64[ns]').astype('float64')
    return vals.flatten().astype('float64')


def _stats(vals):
    valid = vals[~np.isnan(vals)]
    diffs = np.diff(valid)
    s = {
        'n_valid': int(len(valid)),
        'n_nan': int(np.sum(np.isnan(vals))),
    }
    if len(valid):
        s.update(
            min=float(np.min(valid)),
            max=float(np.max(valid)),
            mean=float(np.mean(valid)),
            std=float(np.std(valid)),
        )
    if len(diffs):
        s.update(
            diff_mean=float(np.mean(diffs)),
            diff_std=float(np.std(diffs)),
        )
    return s


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
def test_variable_stats(var):
    vals = _to_float(output[var])
    actual = _stats(vals)
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
def test_variable_diff_stats(var):
    vals = _to_float(output[var])
    actual = _stats(vals)
    exp = expected['variables'][var]
    for stat in ('diff_mean', 'diff_std'):
        if stat not in exp:
            continue
        np.testing.assert_allclose(
            actual[stat], exp[stat], rtol=RTOL, atol=ATOL,
            err_msg=f'{var}: {stat} mismatch',
        )
