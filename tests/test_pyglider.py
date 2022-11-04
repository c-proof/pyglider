import xarray as xr
from pathlib import Path
import pytest
import numpy as np

import pyglider.seaexplorer as seaexplorer
import pyglider.slocum as slocum

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
outname = seaexplorer.raw_to_L0timeseries(rawncdir, l0tsdir,
                                          deploymentyaml, kind='sub')
output = xr.open_dataset(outname)
# Open test data file
test_data = xr.open_dataset(
    library_dir /
    'tests/results/example-seaexplorer/L0-timeseries/dfo-eva035-20190718.nc')
variables = list(output.variables)


def test_variables_seaexplorer():
    test_variables = list(test_data.variables)
    variables.sort()
    test_variables.sort()
    assert variables == test_variables


@pytest.mark.parametrize("var", variables)
def test_example_seaexplorer(var):
    # Test that each variable and its coordinates match
    assert output[var].attrs == test_data[var].attrs
    if var not in ['time']:
        np.testing.assert_allclose(output[var].values, test_data[var].values, rtol=1e-5)
    else:
        dt0 = output[var].values - np.datetime64('2000-01-01')
        dt1 = test_data[var].values - np.datetime64('2000-01-01')
        assert np.allclose(
            np.array(dt0, dtype='float64'),
            np.array(dt1, dtype='float64'))


def test_example_seaexplorer_metadata():
    # Test that attributes match. Have to remove creation and issue dates first
    output.attrs.pop('date_created')
    output.attrs.pop('date_issued')
    test_data.attrs.pop('date_created')
    test_data.attrs.pop('date_issued')
    assert output.attrs == test_data.attrs


# Create an L0 timeseries from slocum data and test that the resulting netcdf is
# identical to the test data
cacdir = str(example_dir / 'example-slocum/cac/') + '/'
sensorlist = str(example_dir / 'example-slocum/dfo-rosie713_sensors.txt')
binarydir = str(example_dir / 'example-slocum/realtime_raw/') + '/'
rawdir_slocum = str(example_dir / 'example-slocum/realtime_rawnc/') + '/'
deploymentyaml_slocum = str(example_dir / 'example-slocum/deploymentRealtime.yml')
l1tsdir = str(example_dir / 'example-slocum/L0-timeseries-test/') + '/'
scisuffix = 'tbd'
glidersuffix = 'sbd'

slocum.binary_to_rawnc(
    binarydir, rawdir_slocum, cacdir, sensorlist, deploymentyaml_slocum,
    incremental=True, scisuffix=scisuffix, glidersuffix=glidersuffix)

slocum.merge_rawnc(rawdir_slocum, rawdir_slocum, deploymentyaml_slocum,
                   scisuffix=scisuffix, glidersuffix=glidersuffix)
outname_slocum = slocum.raw_to_timeseries(
    rawdir_slocum, l1tsdir, deploymentyaml_slocum,
    profile_filt_time=100, profile_min_time=300)
output_slocum = xr.open_dataset(outname_slocum)
# Open test data file
test_data_slocum = xr.open_dataset(
    library_dir /
    'tests/results/example-slocum/L0-timeseries/dfo-rosie713-20190615.nc')
variables_slocum = list(output_slocum.variables)


def test_variables_slocum():
    test_variables = list(test_data_slocum.variables)
    test_variables.sort()
    variables_slocum.sort()
    assert variables_slocum == test_variables


@pytest.mark.parametrize("var", variables_slocum)
def test_example_slocum(var):
    # Test that variables and coordinates match
    assert output_slocum[var].attrs == test_data_slocum[var].attrs
    if var not in ['time']:
        np.testing.assert_allclose(output_slocum[var].values,
                                   test_data_slocum[var].values, rtol=1e-6)
    else:
        dt0 = output_slocum[var].values - np.datetime64('2000-01-01')
        dt1 = test_data_slocum[var].values - np.datetime64('2000-01-01')
        assert np.allclose(
            np.array(dt0, dtype='float64'),
            np.array(dt1, dtype='float64'))


def test_example_slocum_metadata():
    # Test that attributes match. Have to remove creation and issue
    # dates first
    output_slocum.attrs.pop('date_created')
    output_slocum.attrs.pop('date_issued')
    test_data_slocum.attrs.pop('date_created')
    test_data_slocum.attrs.pop('date_issued')
    assert output_slocum.attrs == test_data_slocum.attrs


# Create an L0 timeseries from slocum data and test that the resulting netcdf is
# identical to the test data
cacdir = str(example_dir / 'example-slocum-littleendian/cac/') + '/'
sensorlist = str(
    example_dir / 'example-slocum-littleendian/dfo-maria997_sensors.txt')
binarydir = str(
    example_dir / 'example-slocum-littleendian/realtime_raw/') + '/'
rawdir_slocum = str(
    example_dir / 'example-slocum-littleendian/realtime_rawnc/') + '/'
deploymentyaml_slocum = str(
    example_dir / 'example-slocum-littleendian/deployment.yml')
l1tsdir = str(
    example_dir / 'example-slocum-littleendian/L0-timeseries-test/') + '/'
scisuffix = 'tbd'
glidersuffix = 'sbd'

slocum.binary_to_rawnc(
    binarydir, rawdir_slocum, cacdir, sensorlist, deploymentyaml_slocum,
    incremental=True, scisuffix=scisuffix, glidersuffix=glidersuffix)
slocum.merge_rawnc(rawdir_slocum, rawdir_slocum, deploymentyaml_slocum,
                   scisuffix=scisuffix, glidersuffix=glidersuffix)
outname_slocum_le = slocum.raw_to_timeseries(
    rawdir_slocum, l1tsdir, deploymentyaml_slocum,
    profile_filt_time=400, profile_min_time=100)
output_slocum_le = xr.open_dataset(outname_slocum_le)
# Open test data file
test_data_slocum_le = xr.open_dataset(
    library_dir /
    ('tests/results/example-slocum-littleendian/' +
     'L0-timeseries/dfo-maria997-20220614.nc'))
variables_slocum_le = list(output_slocum.variables)


def test_variables_slocum_littleendian():
    test_variables = list(test_data_slocum_le.variables)
    test_variables.sort()
    variables_slocum.sort()
    assert variables_slocum == test_variables


@pytest.mark.parametrize("var", variables_slocum)
def test_example_slocum_littleendian(var):
    # Test that variables and coordinates match
    assert output_slocum_le[var].attrs == test_data_slocum_le[var].attrs
    if var not in ['time']:
        np.testing.assert_allclose(output_slocum_le[var].values,
                                   test_data_slocum_le[var].values, rtol=1e-6)
    else:
        dt0 = output_slocum_le[var].values - np.datetime64('2000-01-01')
        dt1 = test_data_slocum_le[var].values - np.datetime64('2000-01-01')
        assert np.allclose(
            np.array(dt0, dtype='float64'),
            np.array(dt1, dtype='float64'))


def test_example_slocum_littleendian_metadata():
    # Test that attributes match. Have to remove creation and issue dates first
    output_slocum_le.attrs.pop('date_created')
    output_slocum_le.attrs.pop('date_issued')
    test_data_slocum_le.attrs.pop('date_created')
    test_data_slocum_le.attrs.pop('date_issued')
    assert output_slocum_le.attrs == test_data_slocum_le.attrs
