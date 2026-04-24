# PyGlider: Adjust CTD Variables

PyGlider applies a post-processing protocol to conductivity, temperature, and salinity variables in NetCDF timeseries files, and subsequently generates NetCDF depth–time grids using Python and `xarray`. The resulting NetCDF files are largely CF-compliant.

The basic workflow converts a NetCDF timeseries into an adjusted timeseries and corresponding depth–time grids. This follows the `binary_to_timeseries` (for Slocum gliders) and `raw_to_timeseries` (for Alseamar gliders) protocols, which convert raw glider data into NetCDF format. The variable `outname` refers to the timeseries output from this prior step.

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
  - Or stored in the deployment YAML file
  - If neither is provided, thermal lag correction is not applied

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
ts2 = utils.adjust_CTD(ts, deploymentyaml, interpolate_filter=None)

outfile = f'{ts_path}/{deploy_name}_adjusted.nc'
_log.info('Saving adjusted timeseries to netcdf')
outname = ts2.to_netcdf(outfile)

outname2 = ncprocess.make_gridfiles(
    outname,
    gridpath,
    deploymentyaml,
    fnamesuffix='_adjusted',
    maskfunction=None,
    max_gap=50
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

```python
def flag_CTD_data(
    ts0,
    clean_stdev=3,
    accuracy=None,
):
    """
    Wrapper function to flag CTD data.

    Uses `flag_conductivity_in_depth_space` to flag conductivity as QC1 (good)
    or QC4 (bad) in profile-depth space. Conductivity and salinity are then
    flagged as QC4 wherever conductivity is flagged as QC4.

    Creates `conductivity_QC`, `salinity_QC`, and `temperature_QC` if they do
    not already exist.

    Parameters
    ----------
    ts0 : xarray.Dataset
        Timeseries of mission data.

    clean_stdev : float, optional
        Number of standard deviations from the cleaned mean for data to be
        flagged as QC4.

    Returns
    -------
    ts : xarray.Dataset
        Timeseries of mission data with `conductivity_QC`, `salinity_QC`,
        and `temperature_QC`.
    """
    _log.info('Screening CTD data')

    ts = ts0.copy()

    ts["conductivity"] = ts["conductivity"].where(ts["conductivity"] >= 0.1)

    cond_qc = flag_conductivity_in_depth_space(
            ts,
            d_profile=50,
            dz=5,
            clean_stdev=clean_stdev,
            accuracy=accuracy
    )

    if "conductivity_QC" not in ts.data_vars:
        _log.info('Adding conductivity_QC variable to dataset')

        ts["conductivity_QC"] = xr.DataArray(
            np.ones(ts["conductivity"].shape, dtype=int),
            dims=ts["conductivity"].dims,
            coords=ts["conductivity"].coords,
        )

    if "salinity_QC" not in ts.data_vars:
        _log.info('Adding salinity_QC variable to dataset')
        ts["salinity_QC"] = xr.DataArray(
            np.ones(ts["salinity"].shape, dtype=int),
            dims=ts["salinity"].dims,
            coords=ts["salinity"].coords,
        )

    if "temperature_QC" not in ts.data_vars:
        _log.info('Adding temperature_QC variable to dataset')
        ts["temperature_QC"] = xr.DataArray(
            np.ones(ts["temperature"].shape, dtype=int),
            dims=ts["temperature"].dims,
            coords=ts["temperature"].coords,
        )

    ts["conductivity_QC"] = xr.where(cond_qc == 4, 4, ts["conductivity_QC"])
    ts["salinity_QC"] = xr.where(ts["conductivity_QC"] == 4, 4, ts["salinity_QC"])

    return ts
```

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

```python
def adjust_CTD(
    ts,
    deploymentyaml,
    alpha=None,
    tau=None,
    dTdC=None,
    interpolate_filter=None,
):
    """
    Pulls correction constants from `deploymentyaml`. If `alpha`, `tau`, or `dTdC`
    differ from the values in the YAML file, the values provided as function arguments
    are used and a warning is issued.

    Applies conductivity–temperature lag correction and thermal lag correction when
    the corresponding constants are not `None` or 0. This produces the variables
    `temperature_adjusted` and `salinity_adjusted`. The variables
    `potential_density_adjusted` and `potential_temperature_adjusted` are derived
    from the adjusted temperature and salinity.

    Parameters
    ----------
    ts : xarray.Dataset
        Time series of mission data.

    deploymentyaml : str or list
        Path to a YAML file containing deployment information for the glider.

        If a list is provided, YAML files are read in order, and top-level keys
        in later files overwrite those in earlier files.

    alpha : float, optional
        Thermal lag correction parameter alpha. Default is None.

    tau : float, optional
        Thermal lag correction parameter tau. Default is None.

    dTdC : float, optional
        Time lag (seconds) between temperature and conductivity sensors. Default is None.

    interpolate_filter: callable or None, optional
        Function applied to the dataset before finding the internal temperature.
        Function interpolates over bad data and small data gaps
        to prevent errors from affecting the neighbouring cells. Default is None.
    Returns
    -------
    ts : xarray.Dataset
        Time series dataset with the additional variables:
        `temperature_adjusted`, `salinity_adjusted`,
        `potential_density_adjusted`, and `potential_temperature_adjusted`.
        Metadata are updated to reflect applied corrections.
    """
    logger = logging.getLogger(__name__)
    _log.info('Adjusting CTD data')

    atr = deploymentyaml.get("glider_devices", {}).get("ctd", {})
    thermal = atr.get("Thermal_lag_constants_[alpha,tau]")
    _log.info('CTD thermal lag constants from YAML: %s', thermal)
    yaml_vals = {
        "alpha": thermal[0] if thermal and len(thermal) > 0 else None,
        "tau":   thermal[1] if thermal and len(thermal) > 1 else None,
        "dTdC":  atr.get("dTdC"),
    }

    kw_vals = {
        "alpha": alpha,
        "tau": tau,
        "dTdC": dTdC,
    }

    out = {}
    for key in yaml_vals:
        y = yaml_vals[key]
        k = kw_vals[key]

        if k is not None:
            if y is not None and y != k and logger is not None:
                logger.warning(
                    "%s differs between YAML (%r) and kwargs (%r); using kwargs.",
                    key, y, k
                )
            out[key] = k
        else:
            out[key] = y

    alpha = out.get("alpha", {})
    tau = out.get("tau", {})
    dTdC = out.get("dTdC", {})

    if all(out.get(k) is None for k in ["alpha", "tau", "dTdC"]):
        raise ValueError(
            "Missing required CTD constants after checking kwargs and YAML:'"
            "alpha, tau, dTdC"
        )

    temp_adj = ts.temperature.copy()
    temp_adj.attrs = ts.temperature.attrs.copy()
    temp_adj.attrs["comment"] = "temperature [degC]"

    if dTdC not in (None, 0):
        _log.info('Interpolating temperature data forward by %s seconds', dTdC)
        dt = np.timedelta64(dTdC, "s")
        temp_adj = temp_adj.interp(time=ts.time + dt)

        temp_adj.attrs["history"] = "temperature [degC] adjusted by CT lag"
        temp_adj.attrs["time_lag"] = f"{dTdC} second CT lag corrected"
        ts.attrs["dTdC"] = f"{dTdC} second CT lag corrected"
    else:
        temp_adj.attrs["comment"] = "equivalent to raw temperature"
        ts.attrs["dTdC"] = "No CT lag applied"

    ts["temperature_adjusted"] = temp_adj

    if tau is not None:
        dt = np.diff(ts.time.values).astype("timedelta64[s]").astype(int)
        vals, counts = np.unique(dt, return_counts=True)
        srate = vals[np.argmax(counts)]

        fs = 1 / float(srate)
        fn = 0.5 * fs

        s = apply_thermal_lag(
            ts,
            fn,
            alpha=alpha,
            tau=tau,
            interpolate_filter=interpolate_filter,
        )
        sal_adj = xr.where(ts.salinity_QC == 1, s, ts.salinity)
        sal_adj.attrs = ts.salinity.attrs.copy()
        sal_adj.attrs["history"] = (
            f"adjusted salinity [psu] using thermal lag correction "
            f"(alpha={alpha}, tau={tau})"
        )

        if dTdC not in (None, 0):
            sal_adj.attrs["sources"] = (
                f"conductivity pressure temperature_adjusted "
                f"(corrected for {dTdC} second CT lag)"
            )

        ts["salinity_adjusted"] = sal_adj
        ts.attrs['correction_constants_alpha'] = alpha
        ts.attrs['correction_constants_tau'] = tau

    else:
        _log.info(
            'No thermal lag correction applied; calculating salinity_adjusted '
            'using temperature_adjusted and raw conductivity'
        )
        sal_adj = xr.DataArray(
            gsw.conversions.SP_from_C(
                10 * ts["conductivity"],
                ts["temperature_adjusted"],

                ts.pressure,
            ).values,
            dims=ts.salinity.dims,
            coords=ts.salinity.coords,
        )
        ts.attrs['correction_constants_alpha'] = "None"
        ts.attrs['correction_constants_tau'] = "None"

        sal_adj.attrs = ts.salinity.attrs.copy()

        if dTdC is not None:
            sal_adj.attrs["time_lag"] = "found using temperature_adjusted"

        ts["salinity_adjusted"] = sal_adj

    ts["salinity_adjusted"].attrs = sal_adj.attrs

    ts.attrs["quality_flags"] = (
        "1 = good data; 3 = bad data, potentially correctable; "
        "4 = bad data; 8 = estimated data"
    )

    ts["conductivity"].attrs["comment"] = "raw conductivity"
    ts["conductivity_QC"] = ts.conductivity_QC

    ts["temperature"].attrs["comment"] = "raw temperature [degC]"
    ts["temperature_QC"] = ts.temperature_QC

    ts["temperature_adjusted_QC"] = ts["temperature_QC"]

    ts["salinity"].attrs["comment"] = "raw salinity [psu]"
    ts["salinity_adjusted_QC"] = ts["salinity_QC"]

    ts["density"].attrs["comment"] = "raw density"
    ts["density_QC"] = ts["salinity_QC"]

    ts["potential_density"].attrs["history"] = (
        "calculated using raw salinity and temperature"
    )
    ts["potential_temperature"].attrs["history"] = (
        "calculated using raw salinity and temperature"
    )

    long = ts.longitude.fillna(ts.longitude.mean(skipna=True))
    lat = ts.latitude.fillna(ts.latitude.mean(skipna=True))

    sa_adj = gsw.SA_from_SP(ts["salinity_adjusted"], ts["pressure"], long, lat)
    ct_adj = gsw.CT_from_t(sa_adj, ts["temperature_adjusted"], ts["pressure"])

    _log.info(
        'Calculating potential density and potential temperature using '
        'adjusted salinity and temperature'
    )
    ts["potential_density_adjusted"] = (
        ("time"), 1000 + gsw.density.sigma0(sa_adj, ct_adj).values
    )

    ts["potential_density"].attrs = ts.potential_density.attrs.copy()

    ts["potential_density_adjusted"].attrs["history"] = (
        "calculated using salinity_adjusted and temperature_adjusted"
    )
    ts["potential_density_adjusted"].attrs["sources"] = (
        "salinity_adjusted temperature_adjusted pressure"
    )

    ts["potential_density_adjusted_QC"] = ts["salinity_adjusted_QC"]
    ts["potential_density_adjusted_QC"].attrs = ts.salinity_adjusted_QC.attrs.copy()

    ts["potential_temperature_adjusted"] = (
        ("time"),
        gsw.conversions.pt0_from_t(
            ts.salinity_adjusted,
            ts.temperature_adjusted,
            ts.pressure,
        ).values,
    )

    ts["potential_temperature_adjusted"].attrs = (
        ts.potential_temperature.attrs.copy()
    )
    ts["potential_temperature_adjusted"].attrs["history"] = (
        "calculated using salinity_adjusted and temperature_adjusted"
    )
    ts["potential_temperature_adjusted"].attrs["sources"] = (
        "salinity_adjusted temperature_adjusted pressure"
    )

    ts["potential_temperature_adjusted_QC"] = ts["salinity_adjusted_QC"]
    ts["potential_temperature_adjusted_QC"].attrs = (
        ts.salinity_adjusted_QC.attrs.copy()
    )

    processing_date = date.today().strftime("%Y%m%d")

    vars_ = [
        "salinity_adjusted",
        "temperature_adjusted",
        "potential_density_adjusted",
        "potential_temperature_adjusted",
    ]

    for var in vars_:
        ts[var].attrs["processing_date"] = processing_date

    QC_COMMENT = (
        "1 = good data; 3 = bad data, potentially correctable; "
        "4 = bad data; 8 = estimated data"
    )

    for k in ts.data_vars:
        if k.endswith("_QC"):
            ts[k].attrs["comment"] = QC_COMMENT

    return ts
```

---

### CT lag correction

If `dTdC` is not `None`, temperature is shifted forward in time to align with conductivity.

If `dTdC` is `None` or `0`, no lag correction is applied.

---

### Thermal lag correction

If `alpha` and `tau` are provided, `apply_thermal_lag` is used to:
- Estimate conductivity cell temperature
- Recalculate `salinity_adjusted`

---

## Optional: interpolate_filter

An optional preprocessing step can be applied to reduce noise and prevent spikes from propagating during filtering.

```python
def interpolate_over_salinity_NANs(ds):
    """
    Function applied to the dataset before finding the internal temperature.
    Function interpolates temperature over bad data and small data gaps
    to prevent errors from affecting the neighbouring cells.

    Parameters
    ----------
    ds: DataArray
        Timeseries of mission data

    Returns
    ----------
    interp: DataArray
        Timeseries of interpolated temperature

    """
    _log.info(
        'Interpolating temperature over salinity NaNs and small data gaps '
        'before applying thermal lag correction'
    )
    interp = ds["temperature"].where(ds["temperature_QC"] != 4)
    qc4 = (ds["temperature_QC"] == 4)
    qc4_buf = qc4.rolling(time=5, center=True, min_periods=1).max().astype(bool)
    interp = interp.where(~qc4_buf)

    interp = interp.interpolate_na(
        dim="time",
        method="linear",
        max_gap=np.timedelta64(60, "s"))

    return interp
```

---

## Thermal Lag Function

```python
def apply_thermal_lag(
    ds,
    fn,
    alpha,
    tau,
    interpolate_filter=None,
):
    """
    Function from Garau et al. (2011): estimates temperature inside the
    conductivity cell then recalculates salinity

    Parameters
    ----------
    ds: DataArray
        Timeseries of mission data

    fn: float
        Sampling frequency of the sensor

    alpha : float
        Thermal lag strength constant for the sensor.

    tau: float
        Thermal lag time constant for the sensor.

    interpolate_filter: callable or None, optional
        Function applied to the dataset before finding the internal temperature.
        Function interpolates over bad data and small data gaps
        to prevent errors from affecting the neighbouring cells.

    Returns
    ----------
    sal: DataArray
        Timeseries of salinity_adjusted calculated using the internal temperature of
        the conductivity cell.
    """
    if interpolate_filter is not None:
        temp = interpolate_filter(ds)
        _log.info('Interpolating over bad data and small data gaps before'
                   'applying thermal lag correction')

    else:
        temp = ds.temperature

    _log.info(
        'Applying thermal lag correction with alpha = %s, tau = %s, '
        'and sampling frequency = %s',
        alpha,
        tau,
        fn,
    )
    a = 4 * fn * alpha * tau / (1 + 4*fn*tau)
    b = 1 - 2 * a / alpha
    aa = [1, b]
    bb = [a, -a]
    tempcorr = temp.values.copy()
    tempcell = temp.values.copy()
    good = ~np.isnan(tempcell)
    tempcorr[good] = signal.lfilter(bb, aa, temp.values[good])
    tempcell = tempcell - tempcorr
    sal = gsw.SP_from_C(ds.conductivity * 10, tempcell, ds.pressure)

    return sal
```

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
- QC4 values are excluded from correction steps where appropriate
- Small data gaps can be interpolated prior to filtering to improve stability

---

## Summary

This workflow provides a reproducible method for:
- Flagging bad CTD data
- Correcting sensor response lags
- Producing adjusted physical variables
- Generating gridded NetCDF products

It is designed to be flexible, allowing users to customize correction parameters and preprocessing steps depending on sensor characteristics and mission requirements.

