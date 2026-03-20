# PyGlider: Adjust CTD variables

PyGlider applies a post-processing protocol to conductivity, temperature, and salinity variables in NetCDF timeseries files, and subsequently generates NetCDF depth–time grids using Python and `xarray`. The resulting NetCDF files are largely CF-compliant.

The basic workflow consists of converting a NetCDF timeseries into an adjusted timeseries and corresponding depth–time grids. This follows the `binary_to_timeseries` (for Slocum gliders) and `raw_to_timeseries` (for Alseamar gliders) protocols, which convert raw glider data into NetCDF format. The variable `outname` refers to the timeseries output from this prior step.

```python
outname_ctd = pyglider.ncprocess.adjust_CTD(
    outname,
    deploymentyaml,
    l1tsdir,
    griddir,
    dTdC=None,
    tau=None,
    alpha=None,
    maskfunction=None,
    interp_variables=None
)
```

Data are read from and written to directories, and metadata are supplied via a YAML file.

---

## Post-processing steps within `pyglider.ncprocess.adjust_CTD`

### 1. Identify anomalous conductivity values

We identify and flag conductivity values that are clearly unphysical, typically caused by air bubbles in the conductivity cell.

A two-step statistical criterion is applied:

- First, data points more than **5 standard deviations** from the mean are temporarily flagged within each depth and profile bin.  
- The mean and standard deviation are then recomputed excluding these points.  
- Values still exceeding **3 standard deviations** from the recomputed mean are flagged as **bad (QC = 4)** in `conductivity_QC`.  

However, if the deviation is smaller than the sensor accuracy (0.0003 S/m for the GPCTD), the data are retained.

This procedure is applied using:

- Profile bins of 50 profiles  
- Depth bins of 5 m  

Using profile-index binning (rather than time or temperature) helps isolate unphysical values.

Salinity (`salinity_QC`) is flagged as bad (QC = 4) wherever `conductivity_QC` is QC4.

---

### 2. Determine what dTdC, tau, and alpha are used in the correction

We correct for:

- Sensor misalignment between temperature and conductivity (`dTdC`)  
- Thermal lag effects (`tau`, `alpha`)  

Where:

- `dTdC` = time lag (seconds) between temperature and conductivity sensors  
- `tau` = thermal response time constant (seconds)  
- `alpha` = scaling of thermal coupling between water and the conductivity cell  

Further details and methods for determining these parameters are available at:  
https://cproof.uvic.ca/gliderdata/deployments/reports/

For recent C-PROOF missions, these values are included in the YAML file. However, users may:

- Provide custom values for `dTdC`, `tau`, and `alpha`, or  
- Skip corrections by setting parameters to `None`  

If user-supplied values differ from those in the YAML file, a warning is issued, but the user-provided values are used.

New variables introduced:

- `temperature_adjusted`  
- `salinity_adjusted`  
- `temperature_adjusted_QC`  
- `salinity_adjusted_QC`  

---

### 3. Recalculate derived variables

Using TEOS-10, we recompute:

- `potential_density_adjusted`  
- `potential_temperature_adjusted`  

Their corresponding QC variables:

- `potential_density_adjusted_QC`  
- `potential_temperature_adjusted_QC`  

These are flagged as bad (QC = 4) wherever `salinity_adjusted_QC` is QC4.

---

### 4. Convert the adjusted NetCDF timeseries to a NetCDF depth-time grid

The adjusted timeseries is converted into a gridded NetCDF dataset.

- Default binning:
  - 1 m depth bins  
  - Profile bins  

- Variables are averaged within each bin.

QC variables use a `QC_protocol` that selects the maximum QC value within each bin, ensuring that bad data are not diluted.

Optional parameters:

- `maskfunction`  
- `interp_variables`  

C-PROOF applies:

- `pyglider.ncprocess.CPROOF_mask` to exclude QC4 data from gridded products  
- `pyglider.ncprocess.interpolate_vertical` to interpolate over vertical gaps up to 50 m  
