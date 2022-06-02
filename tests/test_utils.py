import xarray as xr
import numpy as np
import yaml
import pytest
from pathlib import Path
import sys
import pyglider.utils as pgutils

library_dir = Path(__file__).parent.parent.absolute()

test_data = xr.open_dataset(
    library_dir /
    'tests/results/example-seaexplorer/L0-timeseries/dfo-eva035-20190718.nc')
deploymentyaml = (
    library_dir /
    'tests/data-seaexplorer-legato-flntu-arod-ad2cp/deploymentRealtime.yml')
with open(deploymentyaml) as fin:
    deployment = yaml.safe_load(fin)
metadata = deployment['metadata']
ncvar = deployment['netcdf_variables']

ds_unchanged = test_data.copy()
ds_out = pgutils.oxygen_concentration_correction(test_data, ncvar)
vars = list(ds_out)
vars.remove('oxygen_concentration')


@pytest.mark.parametrize("var", vars)
def test_oxygen_concentration_correction_variables_unchanged(var):
    # Test that each variable and its coordinates match
    assert ds_out[var].equals(ds_unchanged[var])


def test_oxygen_concentration_correction():
    oxygen_difference = (np.nanmean(ds_out['oxygen_concentration'].values)
                         - np.nanmean(ds_unchanged['oxygen_concentration'].values))
    assert np.abs(oxygen_difference + 32.328) < 1e-2
