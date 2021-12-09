import xarray as xr
from pathlib import Path
import sys
import pytest

library_dir = Path(__file__).parent.parent.absolute()
sys.path.append(str(library_dir))
import pyglider.seaexplorer as seaexplorer
import pyglider.slocum as slocum

# Create an L0 timeseries from seaexplorer data and test that the resulting netcdf is identical to the test data
rawdir = str(library_dir / 'example-seaexplorer/realtime_raw/') + '/'
rawncdir = str(library_dir / 'example-seaexplorer/realtime_rawnc/') + '/'
deploymentyaml = str(library_dir / 'example-seaexplorer/deploymentRealtime.yml')
l0tsdir = str(library_dir / 'example-seaexplorer/L0-timeseries-test/') + '/'
seaexplorer.raw_to_rawnc(rawdir, rawncdir, deploymentyaml)
seaexplorer.merge_rawnc(rawncdir, rawncdir, deploymentyaml, kind='sub')
outname = seaexplorer.raw_to_L0timeseries(rawncdir, l0tsdir, deploymentyaml, kind='sub')
output = xr.open_dataset(outname)
# Open test data file
test_data = xr.open_dataset(library_dir / 'tests/results/example-seaexplorer/L0-timeseries/dfo-eva035-20190718.nc')
variables = list(output.variables)


def test_variables_seaexplorer():
    test_variables = list(test_data.variables)
    variables.sort()
    test_variables.sort()
    assert variables == test_variables


@pytest.mark.parametrize("var", variables)
def test_example_seaexplorer(var):
    # Test that each variable and its coordinates match
    output_var = output[var].drop('depth')
    test_var = test_data[var].drop('depth')
    assert output_var.equals(test_var)


def test_example_seaexplorer_metadata():
    # Test that attributes match. Have to remove creation and issue dates first
    output.attrs.pop('date_created')
    output.attrs.pop('date_issued')
    test_data.attrs.pop('date_created')
    test_data.attrs.pop('date_issued')
    assert output.attrs == test_data.attrs


# Create an L0 timeseries from slocum data and test that the resulting netcdf is identical to the test data
cacdir = str(library_dir / 'example-slocum/cac/') + '/'
sensorlist = str(library_dir / 'example-slocum/dfo-rosie713_sensors.txt')
binarydir = str(library_dir / 'example-slocum/realtime_raw/') + '/'
rawdir_slocum = str(library_dir / 'example-slocum/realtime_rawnc/') + '/'
deploymentyaml_slocum = str(library_dir / 'example-slocum/deploymentRealtime.yml')
l1tsdir = str(library_dir / 'example-slocum/L0-timeseries-test/') + '/'
scisuffix = 'tbd'
glidersuffix = 'sbd'

slocum.binary_to_rawnc(binarydir, rawdir_slocum, cacdir, sensorlist, deploymentyaml_slocum,
                       incremental=True, scisuffix=scisuffix, glidersuffix=glidersuffix)

slocum.merge_rawnc(rawdir_slocum, rawdir_slocum, deploymentyaml_slocum,
                   scisuffix=scisuffix, glidersuffix=glidersuffix)
outname_slocum = slocum.raw_to_L0timeseries(rawdir_slocum, l1tsdir, deploymentyaml_slocum,
                                            profile_filt_time=100, profile_min_time=300)
output_slocum = xr.open_dataset(outname_slocum)
# Open test data file
test_data_slocum = xr.open_dataset(library_dir / 'tests/results/example-slocum/L0-timeseries/dfo-rosie713-20190615.nc')
variables_slocum = list(output_slocum.variables)


def test_variables_slocum():
    test_variables = list(test_data_slocum.variables)
    test_variables.sort()
    variables_slocum.sort()
    assert variables_slocum == test_variables

@pytest.mark.parametrize("var", variables_slocum)
def test_example_slocum(var):
    # Test that variables and coordinates match
    output_var = output_slocum[var].drop('depth')
    test_var = test_data_slocum[var].drop('depth')
    assert output_var.equals(test_var)


def test_example_slocum_metadata():
    # Test that attributes match. Have to remove creation and issue dates first
    output_slocum.attrs.pop('date_created')
    output_slocum.attrs.pop('date_issued')
    test_data_slocum.attrs.pop('date_created')
    test_data_slocum.attrs.pop('date_issued')
    assert output_slocum.attrs == test_data_slocum.attrs
