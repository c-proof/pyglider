import xarray as xr
from pathlib import Path
import sys

library_dir = Path(__file__).parent.parent.absolute()
sys.path.append(str(library_dir))
import pyglider.seaexplorer as seaexplorer
import pyglider.slocum as slocum


def test_example_seaexplorer():
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
    # Test that variables and coordinates match
    assert output.equals(test_data)
    # Test that attributes match. Have to remove creation and issue dates first
    output.attrs.pop('date_created')
    output.attrs.pop('date_issued')
    test_data.attrs.pop('date_created')
    test_data.attrs.pop('date_issued')
    assert output.attrs == test_data.attrs


def test_example_slocum():
    # Create an L0 timeseries from slocum data and test that the resulting netcdf is identical to the test data
    cacdir = str(library_dir / 'example-slocum/cac/') + '/'
    sensorlist = str(library_dir / 'example-slocum/dfo-rosie713_sensors.txt')
    binarydir = str(library_dir / 'example-slocum/realtime_raw/') + '/'
    rawdir = str(library_dir / 'example-slocum/realtime_rawnc/') + '/'
    deploymentyaml = str(library_dir / 'example-slocum/deploymentRealtime.yml')
    l1tsdir = str(library_dir / 'example-slocum/L0-timeseries-test/') + '/'
    scisuffix = 'tbd'
    glidersuffix = 'sbd'

    slocum.binary_to_rawnc(binarydir, rawdir, cacdir, sensorlist, deploymentyaml,
                           incremental=True, scisuffix=scisuffix, glidersuffix=glidersuffix)

    slocum.merge_rawnc(rawdir, rawdir, deploymentyaml,
                       scisuffix=scisuffix, glidersuffix=glidersuffix)
    outname = slocum.raw_to_L0timeseries(rawdir, l1tsdir, deploymentyaml,
                                         profile_filt_time=100, profile_min_time=300)
    output = xr.open_dataset(outname)
    # Open test data file
    test_data = xr.open_dataset(library_dir / 'tests/results/example-slocum/L0-timeseries/dfo-rosie713-20190615.nc')
    # Test that variables and coordinates match
    assert output.equals(test_data)
    # Test that attributes match. Have to remove creation and issue dates first
    output.attrs.pop('date_created')
    output.attrs.pop('date_issued')
    test_data.attrs.pop('date_created')
    test_data.attrs.pop('date_issued')
    assert output.attrs == test_data.attrs
