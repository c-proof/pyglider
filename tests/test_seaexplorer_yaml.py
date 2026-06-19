"""
YAML-summary-based regression tests for seaexplorer L0 timeseries pipelines.

Intended to replace test_pyglider.py.  Tests cover:
  - variable list, global attributes, per-variable CF attributes
  - per-variable statistics (min, max, mean, std, n_valid, n_nan)
  - per-variable time-derivative statistics (diff_mean, diff_std)
  - time monotonicity
  - NRT sub with interpolation (should not alter values)
  - Raw delayed pipeline with interpolation (must not reduce n_valid)

To regenerate the YAML golden files after an intentional change, run:
    python tests/_generate_expected_yaml.py
"""

import yaml

import numpy as np
import pytest
import xarray as xr
from pathlib import Path

import pyglider.seaexplorer as seaexplorer

library_dir = Path(__file__).parent.parent.absolute()
example_dir = library_dir / 'tests/example-data/'

# ---------------------------------------------------------------------------
# Tolerances and excluded attrs
# ---------------------------------------------------------------------------
from yaml_test_helpers import RTOL, ATOL, SKIP_ATTRS, to_float, stats


# ---------------------------------------------------------------------------
# NRT sub pipeline
# ---------------------------------------------------------------------------
rawdir = str(example_dir / 'example-seaexplorer/realtime_raw/') + '/'
rawncdir = str(example_dir / 'example-seaexplorer/realtime_rawnc/') + '/'
deploymentyaml = str(example_dir / 'example-seaexplorer/deploymentRealtime.yml')
l0tsdir = str(example_dir / 'example-seaexplorer/L0-timeseries-test/') + '/'

seaexplorer.raw_to_rawnc(rawdir, rawncdir, deploymentyaml)
seaexplorer.merge_parquet(rawncdir, rawncdir, deploymentyaml, kind='sub')
outname = seaexplorer.raw_to_L0timeseries(rawncdir, l0tsdir, deploymentyaml, kind='sub')
output = xr.open_dataset(outname)

with open(library_dir / 'tests/expected/example-seaexplorer/L0-timeseries/dfo-eva035-20190718.yml') as f:
    expected = yaml.safe_load(f)


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


def test_time_monotonic():
    t = output['time'].values.astype('datetime64[ns]').astype('float64')
    assert np.all(np.diff(t) > 0), 'time is not monotonically increasing'


# ---------------------------------------------------------------------------
# NRT sub with interpolation — should not alter values for sub data
# ---------------------------------------------------------------------------
import yaml as _yaml  # noqa: E402 (already imported above, just re-alias for clarity)

with open(deploymentyaml) as _fin:
    _dep = _yaml.safe_load(_fin)
_interp_yaml = str(example_dir / 'example-seaexplorer/deploymentRealtimeInterp.yml')
_dep['netcdf_variables']['interpolate'] = True
with open(_interp_yaml, 'w') as _fout:
    _yaml.dump(_dep, _fout)

l0tsdir_interp = str(example_dir / 'example-seaexplorer/L0-timeseries-test-interp/') + '/'
outname_interp = seaexplorer.raw_to_L0timeseries(rawncdir, l0tsdir_interp, _interp_yaml, kind='sub')
output_interp = xr.open_dataset(outname_interp)


@pytest.mark.parametrize('var', sorted(expected['variables']))
def test_interp_nrt_matches_baseline(var):
    """Interpolation on NRT sub data must not change variable statistics."""
    vals = to_float(output_interp[var])
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


# ---------------------------------------------------------------------------
# Raw delayed pipeline
# ---------------------------------------------------------------------------
rawdir_raw = str(example_dir / 'example-seaexplorer-raw/delayed_raw/') + '/'
rawncdir_raw = str(example_dir / 'example-seaexplorer-raw/delayed_rawnc/') + '/'
deploymentyaml_raw = str(example_dir / 'example-seaexplorer-raw/deployment.yml')
l0tsdir_raw = str(example_dir / 'example-seaexplorer-raw/L0-timeseries-test/') + '/'

seaexplorer.raw_to_rawnc(rawdir_raw, rawncdir_raw, deploymentyaml_raw)
seaexplorer.merge_parquet(rawncdir_raw, rawncdir_raw, deploymentyaml_raw, kind='raw')
outname_raw = seaexplorer.raw_to_L0timeseries(
    rawncdir_raw, l0tsdir_raw, deploymentyaml_raw, kind='raw'
)
output_raw = xr.open_dataset(outname_raw)

with open(library_dir / 'tests/expected/example-seaexplorer-raw/L0-timeseries/dfo-bb046-20200908.yml') as f:
    expected_raw = yaml.safe_load(f)


def test_variable_list_raw():
    assert sorted(output_raw.variables) == sorted(expected_raw['variables'])


def test_global_attrs_raw():
    out_attrs = {k: str(v) for k, v in output_raw.attrs.items() if k not in SKIP_ATTRS}
    exp_attrs = {k: v for k, v in expected_raw['attrs'].items() if k not in SKIP_ATTRS}
    assert out_attrs == exp_attrs


@pytest.mark.parametrize('var', sorted(expected_raw['variables']))
def test_variable_attrs_raw(var):
    assert {k: str(v) for k, v in output_raw[var].attrs.items()} == expected_raw['variables'][var]['attrs']


@pytest.mark.parametrize('var', sorted(expected_raw['variables']))
def test_variable_stats_raw(var):
    vals = to_float(output_raw[var])
    actual = stats(vals)
    exp = expected_raw['variables'][var]
    assert actual['n_valid'] == exp['n_valid'], f'{var}: n_valid mismatch'
    assert actual['n_nan'] == exp['n_nan'], f'{var}: n_nan mismatch'
    for stat in ('min', 'max', 'mean', 'std'):
        if stat not in exp:
            continue
        np.testing.assert_allclose(
            actual[stat], exp[stat], rtol=RTOL, atol=ATOL,
            err_msg=f'{var}: {stat} mismatch',
        )


@pytest.mark.parametrize('var', sorted(expected_raw['variables']))
def test_variable_diff_stats_raw(var):
    vals = to_float(output_raw[var])
    actual = stats(vals)
    exp = expected_raw['variables'][var]
    for stat in ('diff_mean', 'diff_std'):
        if stat not in exp:
            continue
        np.testing.assert_allclose(
            actual[stat], exp[stat], rtol=RTOL, atol=ATOL,
            err_msg=f'{var}: {stat} mismatch',
        )


def test_time_monotonic_raw():
    t = output_raw['time'].values.astype('datetime64[ns]').astype('float64')
    assert np.all(np.diff(t) > 0), 'time is not monotonically increasing'


# ---------------------------------------------------------------------------
# Raw delayed with interpolation — must not reduce n_valid counts
# ---------------------------------------------------------------------------
with open(deploymentyaml_raw) as _fin:
    _dep_raw = _yaml.safe_load(_fin)
_interp_yaml_raw = str(example_dir / 'example-seaexplorer-raw/deploymentDelayedInterp.yml')
_dep_raw['netcdf_variables']['interpolate'] = True
with open(_interp_yaml_raw, 'w') as _fout:
    _yaml.dump(_dep_raw, _fout)

l0tsdir_interp_raw = str(example_dir / 'example-seaexplorer-raw/L0-timeseries-test-interp/') + '/'
outname_interp_raw = seaexplorer.raw_to_L0timeseries(
    rawncdir_raw, l0tsdir_interp_raw, _interp_yaml_raw, kind='raw'
)
output_interp_raw = xr.open_dataset(outname_interp_raw)


@pytest.mark.parametrize('var', sorted(expected_raw['variables']))
def test_interp_raw_increases_valid(var):
    """Interpolation on raw delayed data must not reduce non-NaN counts."""
    vals = to_float(output_interp_raw[var])
    n_valid_interp = int(np.sum(~np.isnan(vals)))
    exp_n_valid = expected_raw['variables'][var]['n_valid']
    assert n_valid_interp >= exp_n_valid, (
        f'{var}: interpolated n_valid ({n_valid_interp}) < baseline ({exp_n_valid})'
    )
