# Getting Started: Slocum

## Gather data

Slocum gliders have 4 types of files. For telemetry data there are `*.tbd` files for sensor data, and `*.sbd` for the glider's attitude and position data. These are called `*.ebd` and `*.tbd` respectively, when retrieved from the gliders' payload post deployment. Modern gliders have compressed version of these, eg `*.tcd`, `*.scd` that _pyglider_ should be able to parse. These data files need to be made available in a _single_ directory for _pyglider_ to process. Note that on the glider they are often separated into `science/logs` and `flight/logs`.

Slocum gliders also have a sensor cache file `*.cac`, all of which have randomized names. These are needed by the processing, and are usually stored in a separate cache directory.

You can download example data at <https://cproof.uvic.ca/pyglider-example-data/pyglider-example-data.zip> which will add a local directory `example-data` to your current directory.

## Make a deployment configuration file

The processing routines all take a `deployment.yaml` file as an argument, and information from this is used to fill in metadata and to map sensor names to NetCDF variable names. See {ref}`ExDeplSlocum`, below.

There are four top-levels to the `deployment.yaml`

- `metadata`: The only field that is necessary here is `glider_name`. The rest of the fields will be added to the netcdf files as top-level attributes
- `glider_devices`: This is a list of the glider devices, and any information about them like make, mode, serial number. This is optional, and again is added to the netcdf top-level attributes
- `netcdf_variables`: These are necessary, and map from sensor name (e.g. `source: GPCTD_CONDUCTIVITY`) to a data variable name (e.g. `conductivity`). The fields other than `source:` are optional for the processing to run, and are placed in the attributes of the netCDF variable. However, note that many of these attributes are necessary for CF compliance.
- `profile_variables`: This is a mapping for variables that are per-profile, rather than timeseries. They include variables like a mean position and time for the profile, and a mean derived ocean velocities.

## Process

The example script is relatively straight forward if there is no intermediate processing. See {ref}`ExProcSlocum`, below.

Data comes from an input directory, and is translated into a single CF-compliant netCDF timeseries file using the package [dbdreader](https://dbdreader.readthedocs.io/en/latest/). Finally individual profiles are saved and a 2-D 1-m grid in time-depth is saved.

:::{note}
There is a version that does not require `dbdreader` to do the initial conversion from the Dinkum format to netCDF. However it is quite slow, particularly for full-resolution datasets, and less robust. We suggest using the `slocum.raw_to_timeseries`.
:::

It is possible that between these steps the user will want to add any screening steps, or adjustments to the calibrations. PyGlider does not provide those steps, but is designed so they are easy to add.

(ExDeplSlocum)=

### Example deployment.yaml

```{literalinclude} ../tests/example-data/example-slocum/deploymentRealtime.yml
:language: yaml
```

(ExProcSlocum)=

### Example processing script

```{literalinclude} ../tests/example-data/example-slocum/process_deploymentRealTime.py
:language: python
```
