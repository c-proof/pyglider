# Flexible YAML and OG 1.0 Output

PyGlider 1.0 makes the deployment YAML more flexible so that variable names,
derived-variable computations, and the netCDF dimension name are no longer
hardcoded.  This lets you produce [OceanGliders 1.0 (OG 1.0)](https://oceangliderscommunity.github.io/OG-format-user-manual/OG_Format.html)
trajectory files — with uppercase OG vocabulary names and the `N_MEASUREMENTS`
dimension — using exactly the same processing pipeline as the legacy IOOS GDAC
format.

## Top-level YAML keys

Two new top-level keys control the output format.

### `output_conventions`

Declares the naming convention used in the file.  This value is recorded in
the global attributes and signals to downstream tools which vocabulary the
variable names follow.

```yaml
output_conventions: OG-1.0   # or IOOS_GDAC (default)
```

### `output_dimension`

Sets the name of the observation dimension in the output netCDF file.
Internally pyglider always works with a dimension called `time`; this key
causes it to be renamed on write and restored on read.

```yaml
output_dimension: N_MEASUREMENTS   # default: time
```

If omitted, the dimension name is `time` (IOOS GDAC style).

## Variable flexibility via `processing_role`

In previous versions, the processing pipeline looked for variables by their
literal name — `ds['pressure']`, `ds['conductivity']`, etc. — which forced
the YAML to use those exact keys.

Now each variable entry in `netcdf_variables` can carry a `processing_role`
that tells the pipeline what role the variable plays, regardless of what it is
called in the output file.

```yaml
netcdf_variables:

  PRES:                          # OG 1.0 output name
    source:          sci_water_pressure
    processing_role: pressure    # pipeline looks this up by role
    long_name:       Pressure (measured variable)
    units:           dbar
    ...

  LATITUDE:
    source:          m_lat
    processing_role: latitude
    ...
```

### Required roles

The following roles are used directly by the processing pipeline regardless of
any `processing_method` entries.  If you use non-standard variable names (i.e.
OG 1.0), you **must** declare these roles explicitly:

| Role | Default IOOS GDAC name | Typical OG 1.0 name |
|---|---|---|
| `time` | `time` | `TIME` |
| `latitude` | `latitude` | `LATITUDE` |
| `longitude` | `longitude` | `LONGITUDE` |
| `pressure` | `pressure` | `PRES` |
| `depth` | `depth` | `DEPTH` |
| `profile_index` | `profile_index` | `PROFILE_NUMBER` |

If a required role cannot be resolved to a variable that exists in the dataset,
pyglider raises a `ValueError` with a message pointing to the YAML fix needed.

### Legacy fallback roles

Variables like temperature, conductivity, and oxygen are only looked up by role
in the legacy processing path (when no `processing_method` entries cover
thermodynamic variables).  If you supply `processing_method` blocks for salinity
and density — as OG 1.0 YAMLs should — you do **not** need `processing_role` on
these variables; the method inputs name them explicitly.

If `processing_role` is absent and no `processing_method` covers a variable,
pyglider falls back to looking for a variable whose name matches the role string,
so existing IOOS GDAC YAMLs continue to work without modification.

## Derived variables via `processing_method`

Previously, salinity, density, depth, and profile numbering were computed by
hardcoded calls inside the processing functions, always consuming variables
named `conductivity`, `temperature`, `pressure`, etc.

Now you can specify how each derived variable is computed and which named
inputs to use:

```yaml
  PSAL:
    processing_method:
      practical_salinity:
        conductivity: CNDC
        temperature:  TEMP
        pressure:     PRES
    long_name:    Sea water practical salinity
    units:        "1"
    ...

  DEPTH:
    processing_method:
      depth_from_pressure:
        pressure: PRES
        latitude: LATITUDE
    processing_role: depth
    ...

  PROFILE_NUMBER:
    processing_method:
      find_profiles:
        pressure: PRES
    processing_role: profile_index
    ...
```

The `processing_method` key contains a single-entry mapping from a method name
to its named inputs.  The inputs are references to other variable names in
`netcdf_variables`.

### Built-in method names

| Method | Computes | Required inputs |
|---|---|---|
| `practical_salinity` | SP via TEOS-10 | `conductivity`, `temperature`, `pressure` |
| `potential_temperature` | θ via TEOS-10 | `salinity`, `temperature`, `pressure` |
| `potential_density_sigma0` | σ₀ via TEOS-10 | `salinity`, `temperature`, `pressure`, `latitude`, `longitude` |
| `density` | in-situ density | `salinity`, `temperature`, `pressure`, `latitude`, `longitude` |
| `depth_from_pressure` | depth (m) | `pressure`, `latitude` |
| `find_profiles` | profile index and direction | `pressure` |
| `distance_over_ground` | cumulative distance | `latitude`, `longitude` |

### Custom methods

If the method name contains a `.` it is treated as a dotted Python import
path.  The function must have the signature:

```python
def my_method(ds: xr.Dataset, inputs: dict, output_name: str) -> xr.DataArray:
    ...
```

where `inputs` is `{role_name: variable_name_in_ds, ...}` as declared in the
YAML.

## Dimension handling

The `time` dimension is always used internally.  On write, if
`output_dimension` is set to something other than `time`, the dimension is
renamed just before the file is saved.  On read (e.g. when loading a
timeseries file to make grids or profiles), pyglider detects the non-standard
dimension by finding the coordinate with `standard_name: time` and renames it
back to `time`.  This is transparent to any intermediate processing steps.

## OG 1.0 example

See the complete example YAML and processing script for an OG 1.0 Slocum
deployment:

```{literalinclude} ../tests/example-data/example-slocum/deploymentRealtime_og10.yml
:language: yaml
```

```{literalinclude} ../tests/example-data/example-slocum/process_deploymentRealtime_og10.py
:language: python
```
