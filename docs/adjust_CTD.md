# PyGlider: Adjust CTD Variables

PyGlider applies a post-processing protocol to conductivity, temperature, and salinity variables in NetCDF timeseries files, and then generates depth–time NetCDF grids using `xarray`. The resulting NetCDF files are largely CF-compliant.

The basic workflow converts a NetCDF timeseries into an adjusted timeseries and corresponding depth–time grids. This follows the `binary_to_timeseries` (for Slocum gliders) and `raw_to_timeseries` (for Alseamar gliders) protocols, which convert raw glider data into NetCDF format. outname refers to the NetCDF timeseries produced by the initial conversion step. 

The adjusted file is generated through three steps: flagging CTD data, adjusting CTD variables, and gridding the time series. This workflow uses known thermal lag constants and the lag between temperature and conductivity signals for each sensor.

An example of how to determine these constants is available at:
https://cproof.uvic.ca/gliderdata/deployments/reports/

---

## Workflow Overview

The CTD adjustment pipeline consists of three main steps:

1. **Flag CTD data**  
   Identify and flag unphysical conductivity, temperature, and salinity values.

2. **Apply CTD corrections**  
   - Correct temperature–conductivity lag (`dTdC`)
   - Apply thermal lag correction (`alpha`, `tau`)

3. **Generate gridded products**  
   Create depth–time NetCDF grids from the adjusted time series.

---

## User Options

The following parameters can be customized:

- **alpha, tau**  
  Thermal lag correction constants.
  - Can be provided as function arguments
  - Or stored in the deployment YAML file or an new YML file. Note that `utils._get_deployment(deploymentyaml)` can take a list of files in *deploymentyaml* and parse them
    for deployment information, with subsequent files overwriting previous files.
  - If information is provided in either method, a thermal lag correction is not applied.

- **dTdC**  
  Time lag (seconds) between temperature and conductivity sensors

- **interpolate_filter**  
  Optional function to interpolate over small gaps before applying thermal lag correction

- **max_gap (in gridding)**  
  Maximum vertical gap size (in meters) to interpolate

---

## Example Processing Script

```python
import logging
from pathlib import Path
import xarray as xr
import pyglider.ncprocess as ncprocess
import pyglider.utils as utils

logging.basicConfig(level=logging.INFO)
_log = logging.getLogger(__name__)

cwd = Path.cwd()
deploy_name = cwd.name
glider_name = cwd.parent.name

openfile = f'./L0-timeseries/{deploy_name}_delayed.nc'
deploymentyaml = './deploymentRealtime.yml'
gridpath = './L0-gridfiles/'
ts_path = './L0-timeseries/'

ts = xr.open_dataset(openfile)
deployment = utils._get_deployment(deploymentyaml)
ts = utils.flag_CTD_data(ts)
ts2 = utils.adjust_CTD(ts, deployment, interpolate_filter=None)

outfile = f'{ts_path}/{deploy_name}_adjusted.nc'
_log.info('Saving adjusted timeseries to netcdf')
outname = ts2.to_netcdf(outfile)

outname2 = ncprocess.make_gridfiles(
    outname,
    gridpath,
    deploymentyaml,
    fnamesuffix='_adjusted',
    maskfunction=None,
    max_gap=100
)
```

The following sections describe each processing step in detail.

---

## flag_CTD_data

This step identifies and flags CTD values that are clearly unphysical (QC4), typically caused by air bubbles in the conductivity cell.

Conductivity data are grouped into profile bins (`d_profile`) and depth bins (`dz`). Within each depth bin:

1. Values more than 5 standard deviations from the mean are temporarily excluded
2. A cleaned mean and standard deviation are recomputed
3. Values exceeding `clean_stdev` from the cleaned mean are flagged as QC4

If `accuracy` is provided, small deviations below this threshold are not flagged.

General behavior:
- QC variables are added if missing
- Values not flagged as QC4 are set to QC1
- `salinity_QC` is set to QC4 wherever `conductivity_QC` is QC4

---

## adjust_CTD

This step applies two corrections:

1. **CT lag correction (`dTdC`)**
2. **Thermal lag correction (`alpha`, `tau`)**

Correction constants can be:
- Passed as function arguments
- Loaded from the deployment YAML file

If both are provided and differ, function arguments take precedence and a warning is issued.

### Example YAML configuration

```yaml
glider_devices:
  ctd:
    Thermal_lag_constants_[alpha,tau]: [0.2, 2]
    dTdC: 0
```

---

### CT lag correction

If `dTdC` is not `None`, temperature is shifted back in time to align with conductivity.

If `dTdC` is `None` or `0`, no lag correction is applied.

---

### Thermal lag correction

If `alpha` and `tau` are provided, `apply_thermal_lag` is used to:
- Estimate conductivity cell temperature
- Recalculate `salinity_adjusted`

---

## Optional: interpolate_filter

An optional preprocessing step can be applied to reduce noise and prevent spikes from propagating when applying the thermal lag correction. This function should take and return an xarray Dataset of salinity (e.g., linear interpolation over short gaps). 

Example: `utils.interpolate_over_salinity_NANs` 

---

## Output Variables

The adjusted dataset includes:

- `temperature_adjusted`
- `salinity_adjusted`
- `potential_density_adjusted`
- `potential_temperature_adjusted`

Quality control flags are propagated to adjusted variables.

Metadata are updated to document:
- Applied corrections
- Processing date
- Data provenance

---

## Notes

- Thermal lag correction is only applied if `alpha` and `tau` are defined
- CT lag correction is only applied if `dTdC` is non-zero
- QC4 values may be excluded from processing if a masking function is applied during gridding
- Small data gaps can be interpolated prior to filtering to improve stability

---

## Summary

This workflow provides a reproducible method for:
- Flagging bad CTD data
- Correcting sensor response lags
- Producing adjusted physical variables
- Generating gridded NetCDF products

It is designed to be flexible, allowing users to customize correction parameters and preprocessing steps depending on sensor characteristics and mission requirements.

