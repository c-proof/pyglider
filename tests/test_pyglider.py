from pathlib import Path

import numpy as np
import pytest
import xarray as xr
import yaml

import pyglider.seaexplorer as seaexplorer

library_dir = Path(__file__).parent.parent.absolute()
example_dir = library_dir / 'tests/example-data/'

# Create an L0 timeseries from seaexplorer data and test that the resulting netcdf
# is identical to the test data
rawdir = str(example_dir / 'example-seaexplorer/realtime_raw/') + '/'
rawncdir = str(example_dir / 'example-seaexplorer/realtime_rawnc/') + '/'
deploymentyaml = str(example_dir / 'example-seaexplorer/deploymentRealtime.yml')
l0tsdir = str(example_dir / 'example-seaexplorer/L0-timeseries-test/') + '/'
seaexplorer.raw_to_rawnc(rawdir, rawncdir, deploymentyaml)
seaexplorer.merge_parquet(rawncdir, rawncdir, deploymentyaml, kind='sub')
outname = seaexplorer.raw_to_L0timeseries(rawncdir, l0tsdir, deploymentyaml, kind='sub')
output = xr.open_dataset(outname)
# Open test data file
test_data = xr.open_dataset(
    library_dir
    / 'tests/expected/example-seaexplorer/L0-timeseries/dfo-eva035-20190718.nc'
)
variables = list(output.variables)


def test_variables_seaexplorer():
    test_variables = list(test_data.variables)
    variables.sort()
    test_variables.sort()
    assert variables == test_variables


@pytest.mark.parametrize('var', variables)
def test_example_seaexplorer(var):
    # Test that each variable and its coordinates match
    assert output[var].attrs == test_data[var].attrs
    if var not in ['time']:
        np.testing.assert_allclose(output[var].values, test_data[var].values, rtol=1e-5)
    else:
        dt0 = output[var].values - np.datetime64('2000-01-01')
        dt1 = test_data[var].values - np.datetime64('2000-01-01')
        assert np.allclose(
            np.array(dt0, dtype='float64'), np.array(dt1, dtype='float64')
        )


def test_example_seaexplorer_metadata():
    # Test that attributes match. Have to remove creation and issue dates first
    output.attrs.pop('date_created')
    output.attrs.pop('date_issued')
    test_data.attrs.pop('date_created')
    test_data.attrs.pop('date_issued')
    assert output.attrs == test_data.attrs


# Test that interpolation over nans does not change the output with nrt data
with open(deploymentyaml) as fin:
    deployment = yaml.safe_load(fin)
interp_yaml = str(example_dir / 'example-seaexplorer/deploymentRealtimeInterp.yml')
deployment['netcdf_variables']['interpolate'] = True
with open(interp_yaml, 'w') as fout:
    yaml.dump(deployment, fout)
l0tsdir_interp = (
    str(example_dir / 'example-seaexplorer/L0-timeseries-test-interp/') + '/'
)

outname_interp = seaexplorer.raw_to_L0timeseries(
    rawncdir, l0tsdir_interp, interp_yaml, kind='sub'
)
output_interp = xr.open_dataset(outname_interp)


@pytest.mark.parametrize('var', variables)
def test_example_seaexplorer_interp_nrt(var):
    assert output_interp[var].attrs == test_data[var].attrs
    if var not in ['time']:
        np.testing.assert_allclose(
            output_interp[var].values, test_data[var].values, rtol=1e-5
        )
    else:
        dt0 = output_interp[var].values - np.datetime64('2000-01-01')
        dt1 = test_data[var].values - np.datetime64('2000-01-01')
        assert np.allclose(
            np.array(dt0, dtype='float64'), np.array(dt1, dtype='float64')
        )


# Test raw (full resolution) seaexplorer data.
rawdir = str(example_dir / 'example-seaexplorer-raw/delayed_raw/') + '/'
rawncdir = str(example_dir / 'example-seaexplorer-raw/delayed_rawnc/') + '/'
deploymentyaml_raw = str(example_dir / 'example-seaexplorer-raw/deployment.yml')
l0tsdir = str(example_dir / 'example-seaexplorer-raw/L0-timeseries-test/') + '/'
seaexplorer.raw_to_rawnc(rawdir, rawncdir, deploymentyaml_raw)
seaexplorer.merge_parquet(rawncdir, rawncdir, deploymentyaml_raw, kind='raw')
outname_raw = seaexplorer.raw_to_L0timeseries(
    rawncdir, l0tsdir, deploymentyaml_raw, kind='raw', deadreckon=True
)
output_raw = xr.open_dataset(outname_raw)
# Open test data file
test_data_raw = xr.open_dataset(
    library_dir
    / 'tests/expected/example-seaexplorer-raw/L0-timeseries/dfo-bb046-20200908.nc'
)


@pytest.mark.parametrize('var', variables)
def test_example_seaexplorer_raw(var):
    # Test that each variable and its coordinates match
    assert output_raw[var].attrs == test_data_raw[var].attrs
    if var not in ['time']:
        np.testing.assert_allclose(
            output_raw[var].values, test_data_raw[var].values, rtol=1e-5
        )
    else:
        dt0 = output_raw[var].values - np.datetime64('2000-01-01')
        dt1 = test_data_raw[var].values - np.datetime64('2000-01-01')
        assert np.allclose(
            np.array(dt0, dtype='float64'), np.array(dt1, dtype='float64')
        )


def test_example_seaexplorer_metadata_raw():
    # Test that attributes match. Have to remove creation and issue dates first
    output_raw.attrs.pop('date_created')
    output_raw.attrs.pop('date_issued')
    test_data_raw.attrs.pop('date_created')
    test_data_raw.attrs.pop('date_issued')
    assert output_raw.attrs == test_data_raw.attrs


# Test nan interpolation on raw data.

# Test that interpolation over nans in raw data results in a greater or equal number of non-nan values
with open(deploymentyaml_raw) as fin:
    deployment_raw = yaml.safe_load(fin)
interp_yaml = str(example_dir / 'example-seaexplorer-raw/deploymentDelayedInterp.yml')
deployment_raw['netcdf_variables']['interpolate'] = True
with open(interp_yaml, 'w') as fout:
    yaml.dump(deployment_raw, fout)
l0tsdir_interp_raw = (
    str(example_dir / 'example-seaexplorer-raw/L0-timeseries-test-interp/') + '/'
)

outname_interp_raw = seaexplorer.raw_to_L0timeseries(
    rawncdir, l0tsdir_interp_raw, interp_yaml, kind='raw', deadreckon=True
)
output_interp_raw = xr.open_dataset(outname_interp_raw)


@pytest.mark.parametrize('var', variables)
def test_example_seaexplorer_interp_raw(var):
    assert output_interp_raw[var].attrs == test_data_raw[var].attrs
    if var not in ['time']:
        assert np.count_nonzero(
            ~np.isnan(output_interp_raw[var].values)
        ) >= np.count_nonzero(~np.isnan(test_data_raw[var].values))
    else:
        dt0 = output_interp_raw[var].values - np.datetime64('2000-01-01')
        dt1 = test_data_raw[var].values - np.datetime64('2000-01-01')
        assert np.allclose(
            np.array(dt0, dtype='float64'), np.array(dt1, dtype='float64')
        )
