from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from pyglider.process_adjusted import run_process_adjusted

library_dir = Path(__file__).parent.parent.absolute()
example_dir = library_dir / 'tests/expected/example-seaexplorer'
expected_dir = library_dir / 'tests/expected/example-seaexplorer'
yaml_dir = library_dir / 'tests/example-data/example-seaexplorer'

outname = run_process_adjusted(
    example_dir,
    deploy_name = 'dfo-eva035-20190718',
    deployfile=yaml_dir / 'deploymentRealtime.yml',
    adjustedyaml=yaml_dir / 'adjusted.yml',
)
output = xr.open_dataset(outname)

test_data = xr.open_dataset(
    expected_dir / 'L0-timeseries' / 'dfo-eva035-20190718_adjusted.nc'
)

variables = list(output.variables)


def test_variables_process_adjusted():
    test_variables = list(test_data.variables)
    variables.sort()
    test_variables.sort()
    assert variables == test_variables


@pytest.mark.parametrize('var', variables)
def test_process_adjusted_timeseries(var):
    assert output[var].attrs == test_data[var].attrs
    if var not in ['time']:
        np.testing.assert_allclose(
            output[var].values,
            test_data[var].values,
            rtol=1e-5,
            equal_nan=True,
        )
    else:
        dt0 = output[var].values - np.datetime64('2000-01-01')
        dt1 = test_data[var].values - np.datetime64('2000-01-01')
        assert np.allclose(
            np.array(dt0, dtype='float64'),
            np.array(dt1, dtype='float64'),
        )


def test_process_adjusted_metadata():
    output.attrs.pop('date_created', None)
    output.attrs.pop('date_issued', None)
    test_data.attrs.pop('date_created', None)
    test_data.attrs.pop('date_issued', None)
    assert output.attrs == test_data.attrs
