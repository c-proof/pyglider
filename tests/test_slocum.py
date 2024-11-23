import xarray as xr

from compliance_checker.runner import ComplianceChecker, CheckSuite
import json
from pathlib import Path
import pytest
import numpy as np
import yaml

import pyglider.ncprocess as ncprocess
import pyglider.seaexplorer as seaexplorer
import pyglider.slocum as slocum



library_dir = Path(__file__).parent.parent.absolute()
example_dir = library_dir / 'tests/example-data/'

# Create an L0 timeseries from slocum data and test that the resulting netcdf is
# identical to the test data
cacdir = example_dir / 'example-slocum/cac/'
sensorlist = str(example_dir / 'example-slocum/dfo-rosie713_sensors.txt')
binarydir = str(example_dir / 'example-slocum/realtime_raw/') + '/'
rawdir_slocum = str(example_dir / 'example-slocum/realtime_rawnc/') + '/'
deploymentyaml_slocum = str(example_dir / 'example-slocum/deploymentRealtime.yml')
tsdir = str(example_dir / 'example-slocum/L0-timeseries/') + '/'
scisuffix = 'tbd'
glidersuffix = 'sbd'
profiledir = str(example_dir /  'example-slocum/L0-profiles/')
do_direct = True

# This needs to get run every time the tests are run, so do at top level:

# turn *.sbd and *.tbd into timeseries netcdf files
outname_slocum = slocum.binary_to_timeseries(binarydir, cacdir, tsdir, deploymentyaml_slocum,
                                             search='*.[s|t]bd', profile_filt_time=20,
                                             profile_min_time=20)
# make profiles...
ncprocess.extract_timeseries_profiles(outname_slocum, profiledir, deploymentyaml_slocum,
                                      force=True)

output_slocum = xr.open_dataset(outname_slocum)
# Open test data file
test_data_slocum = xr.open_dataset(
    library_dir /
    'tests/expected/example-slocum/L0-timeseries/dfo-rosie713-20190615.nc')
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
    output_slocum.attrs.pop('history')
    output_slocum.attrs.pop('netcdf_c_version')
    output_slocum.attrs.pop('netcdf_py_version')
    test_data_slocum.attrs.pop('date_created')
    test_data_slocum.attrs.pop('date_issued')
    test_data_slocum.attrs.pop('history')
    test_data_slocum.attrs.pop('netcdf_version')
    assert output_slocum.attrs == test_data_slocum.attrs


# test the profiles with compliance_checker...

def test_profiles_compliant():
    # Load all available checker classes
    check_suite = CheckSuite()
    check_suite.load_all_available_checkers()
    # Run cf and adcc checks
    path = profiledir + '/dfo-rosie713-20190620T1313.nc'
    checker_names = ['gliderdac', 'cf:1.8']
    verbose = 0
    criteria = 'normal'
    output_filename = example_dir / 'report.json'
    output_format = 'json'
    """
    Inputs to ComplianceChecker.run_checker

    path            Dataset location (url or file)
    checker_names   List of string names to run,
                    should match keys of checkers dict (empty list means run all)
    verbose         Verbosity of the output (0, 1, 2)
    criteria        Determines failure (lenient, normal, strict)
    output_filename Path to the file for output
    output_format   Format of the output

    @returns                If the tests failed (based on the criteria)
    """
    return_value, errors = ComplianceChecker.run_checker(path,
                                                        checker_names,
                                                        verbose,
                                                        criteria,
                                                        output_filename=output_filename,
                                                        output_format=output_format)
    # Open the JSON output and get the compliance scores
    with open(output_filename, 'r') as fp:
        cc_data = json.load(fp)
        test = cc_data['gliderdac']
        assert test['high_count'] == 0
        assert test['medium_count'] == 0
        assert test['low_count'] == 0
        test = cc_data['cf:1.8']
        assert test['high_count'] == 0
        assert test['medium_count'] == 0
        assert test['low_count'] == 0


def test_timeseries_compliant():
    # Load all available checker classes
    check_suite = CheckSuite()
    check_suite.load_all_available_checkers()
    # Run cf and adcc checks
    path = tsdir + '/dfo-rosie713-20190615.nc'
    checker_names = ['cf:1.8']
    verbose = 0
    criteria = 'normal'
    output_filename = example_dir / 'report.json'
    output_format = 'json'
    """
    Inputs to ComplianceChecker.run_checker

    path            Dataset location (url or file)
    checker_names   List of string names to run,
                    should match keys of checkers dict (empty list means run all)
    verbose         Verbosity of the output (0, 1, 2)
    criteria        Determines failure (lenient, normal, strict)
    output_filename Path to the file for output
    output_format   Format of the output

    @returns                If the tests failed (based on the criteria)
    """
    return_value, errors = ComplianceChecker.run_checker(path,
                                                        checker_names,
                                                        verbose,
                                                        criteria,
                                                        output_filename=output_filename,
                                                        output_format=output_format)
    # Open the JSON output and get the compliance scores
    with open(output_filename, 'r') as fp:
        cc_data = json.load(fp)
        test = cc_data['cf:1.8']
        assert test['high_count'] == 0
        # somehow the checker is confused by our trajectory variables.
        assert test['medium_count'] == 1
        assert test['low_count'] == 0
