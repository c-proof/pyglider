# Getting Started: SeaExplorer


## Gather data

SeaExplorers send back and record two main types of files, glider files (`*.gli.*`) that contain glider navigation information, and payload files (`*.pld1.*`) that contain the science data.  These can be subset files, `*.sub.*` that Alseamar decimates for transmission, or they can be full resolution files from the glider (`*.raw.*`), offloaded post mission.  The raw or subset files need to be made available in a single directory for `pyglider` to process.

You can download and expand example data using `.get_example_data`:

```python
import pyglider.example_data as pexamp

pexamp.get_example_data('./')
```

which will add a local directory `example-data` to your current directory.

## Make a deployment configuration file

The processing routines all take a `deployment.yaml` file as an argument, and information from this is used to fill in metadata and to map sensor names to NetCDF variable names.  See {ref}`ExDepl`, below.

There are four top-levels to the `deployment.yaml`

- `metadata`: The only field that is necessary here is `glider_name`.  The rest of the fields will be added to the netcdf files as top-level attributes
- `glider_devices`: This is a list of the glider devices, and any information about them like make, mode, serial number.  This is optional, and again is added to the netcdf top-level attributes
- `netcdf_variables`: These are necessary, and map from sensor name (e.g. `source: GPCTD_CONDUCTIVITY`) to a data variable name (e.g. `conductivity`).  The fields other than `source:` are optional for the processing to run, and are placed in the attributes of the netCDF variable.  However, note that many of these attributes are necessary for CF compliance.
- `profile_variables`: This is a mapping for variables that are per-profile, rather than timeseries.  They include variables like a mean position and time for the profile, and a mean derived ocean velocities.

## Process

The example script is relatively straight forward if there is no intermediate processing.  See {ref}`ExProc`, below.

Data comes from an input directory, and is translated to raw glider-dependent parquet files files and put in a new directory. These files are useful of their own right. Apache Parquet is a columnar oriented format for storing tabular data. Parquet files take up less space than netCDF or csv and are much faster to read and write. These files can be opened with [polars.read_parquet](https://pola-rs.github.io/polars-book/user-guide/howcani/io/parquet.html) or [pandas.read_parquet](https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.read_parquet.html). These files are then merged into a single monolithic parquet file, and this is translated to a CF-compliant timeseries netCDF file.  Finally individual profiles are saved and a 2-D 1-m grid in time-depth is saved.

It is likely that between these steps the user will want to add any screening steps, or adjustments to the calibrations.  PyGlider does not provide those steps.


(ExDepl)=
### Example deployment.yaml

```{literalinclude} ../tests/example-data/example-seaexplorer/deploymentRealtime.yml
:language: yaml
```

(ExProc)=
### Example processing script

```{literalinclude}  ../tests/example-data/example-seaexplorer/process_deploymentRealTime.py
:language: python
```
