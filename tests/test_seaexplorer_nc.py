"""
CDL/NC-based regression tests for seaexplorer L0 timeseries pipelines.

Golden files live in tests/expected/ as cleaned .nc files (SKIP_ATTRS stripped)
with .cdl headers alongside for git readability.  Tests compare pipeline output
against golden NC using xarray.testing.assert_identical.

To regenerate golden files after an intentional change, run:
    python tests/_generate_expected_cdl.py
"""

from pathlib import Path

import numpy as np
import pytest
import xarray as xr
import yaml

import pyglider.seaexplorer as seaexplorer

from nc_test_helpers import assert_datasets_equal

library_dir = Path(__file__).parent.parent.absolute()
example_dir = library_dir / 'tests/example-data/'
expected_dir = library_dir / 'tests/expected/'

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
output = xr.open_dataset(outname).load()

expected_nc = xr.open_dataset(
    expected_dir / 'example-seaexplorer/L0-timeseries/dfo-eva035-20190718.nc'
).load()


def test_dataset():
    assert_datasets_equal(output, expected_nc)


def test_time_monotonic():
    t = output['time'].values.astype('datetime64[ns]').astype('float64')
    assert np.all(np.diff(t) > 0), 'time is not monotonically increasing'


# ---------------------------------------------------------------------------
# NRT sub with interpolation — should not alter values for sub data
# ---------------------------------------------------------------------------
with open(deploymentyaml) as _fin:
    _dep = yaml.safe_load(_fin)
_interp_yaml = str(example_dir / 'example-seaexplorer/deploymentRealtimeInterp.yml')
_dep['netcdf_variables']['interpolate'] = True
with open(_interp_yaml, 'w') as _fout:
    yaml.dump(_dep, _fout)

l0tsdir_interp = str(example_dir / 'example-seaexplorer/L0-timeseries-test-interp/') + '/'
outname_interp = seaexplorer.raw_to_L0timeseries(rawncdir, l0tsdir_interp, _interp_yaml, kind='sub')
output_interp = xr.open_dataset(outname_interp).load()

def test_interp_nrt_matches_baseline():
    """Interpolation on NRT sub data must not change data values."""
    xr.testing.assert_equal(output_interp, expected_nc)


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
output_raw = xr.open_dataset(outname_raw).load()

expected_raw_nc = xr.open_dataset(
    expected_dir / 'example-seaexplorer-raw/L0-timeseries/dfo-bb046-20200908.nc'
).load()


def test_dataset_raw():
    assert_datasets_equal(output_raw, expected_raw_nc)


def test_time_monotonic_raw():
    t = output_raw['time'].values.astype('datetime64[ns]').astype('float64')
    assert np.all(np.diff(t) > 0), 'time is not monotonically increasing'


# ---------------------------------------------------------------------------
# Raw delayed with interpolation — must not reduce n_valid counts
# ---------------------------------------------------------------------------
with open(deploymentyaml_raw) as _fin:
    _dep_raw = yaml.safe_load(_fin)
_interp_yaml_raw = str(example_dir / 'example-seaexplorer-raw/deploymentDelayedInterp.yml')
_dep_raw['netcdf_variables']['interpolate'] = True
with open(_interp_yaml_raw, 'w') as _fout:
    yaml.dump(_dep_raw, _fout)

l0tsdir_interp_raw = str(example_dir / 'example-seaexplorer-raw/L0-timeseries-test-interp/') + '/'
outname_interp_raw = seaexplorer.raw_to_L0timeseries(
    rawncdir_raw, l0tsdir_interp_raw, _interp_yaml_raw, kind='raw'
)
output_interp_raw = xr.open_dataset(outname_interp_raw).load()


def _n_valid(da: xr.DataArray) -> int:
    vals = da.values
    if not np.issubdtype(vals.dtype, np.number) and not np.issubdtype(vals.dtype, np.datetime64):
        return 0
    if np.issubdtype(vals.dtype, np.datetime64):
        vals = vals.astype('float64')
    return int(np.sum(~np.isnan(vals.flatten().astype('float64'))))


@pytest.mark.parametrize('var', sorted(expected_raw_nc.variables))
def test_interp_raw_increases_valid(var):
    """Interpolation on raw delayed data must not reduce non-NaN counts."""
    n_interp = _n_valid(output_interp_raw[var])
    n_base = _n_valid(expected_raw_nc[var])
    assert n_interp >= n_base, (
        f'{var}: interpolated n_valid ({n_interp}) < baseline ({n_base})'
    )
