
# Slocum processing steps

## Setup:

- install `pyglider` and make sure it can be imported.
- remember to activate your conda environment if you are using
  anaconda to keep environments straight.  i.e. I do `conda activate pyglider`

## Get raw data:

For these examples its being put in `realtime_raw`, though you can change
that name, and should for post-processing data.  In this case, we are getting
the data from Teledyne Webb directly by syncing our data directories with
our computer.  The raw data could also come from a card offload.

The files have different suffixes if they come from the realtime data:
`*.sbd` are the glider data and `*.tbd` are the science data.  For
post-processing these would be `dbd` and `edb` respectively.

All the files we want should go in here.  If you dont want to process all the
files you can manually remove them.  

### Get the headers:

The slocums need to know the layout of the data in the files from `*.cac`
files.  These often change if the glider changes setups, and need to be
offloaded from the glider.  Put these files in `./cac/` so the processing
can find the right files.

Or get these from SFMC: `Configuration/Cache Files` and put them in `cac`


### Get and edit the sensor filter file:

```
./dfo_rosie713_sensors.txt
```
This is the list of data to try and parse out of the slocum files.  These will
get translated to netcdf files at the next steps.


## Edit deployment info:

Edit this to have correct glider name and as much meta info as needed.  
`deploymentRealtime.yml`

## Edit `process_deploymentRealtime.py`

Note that you will want to edit this to have correct directories:
Also, realtime and delayed time have different suffixies:

```
scisuffix    = 'tbd'
glidersuffix = 'sbd'
```

## Run `python process_deploymentRealtime.py`

This should do a number of steps:

1. convert the `*.tdb/sbd` files to raw netcdf files (`realtime_rawnc/`) and then concatenate them to one or more merged files.  For instance in the
example data, there are four files created: `dfo-rosie713-0097-rawdbd.nc`,
`dfo-rosie713-0097-rawebd.nc`, `dfo-rosie713-0098-rawdbd.nc` and `dfo-rosie713-0098-rawebd.nc`, that are a concatenation of the missions.  These are not
concatenated into one file because different missions sometimes have
different sensors activated and they won't merge smoothly.  

2. The science files are created in `L1_timeseries`, `L1_profiles`, and
`L2_gridfiles`. These are proper nectdf files that have metadata and attributes, and should be compliant with US-IOOS standards.

3. A plot will be made in `figs/`.  
