
PyGlider: Convert glider data to NetCDF
=======================================

PyGlider converts raw glider files to NetCDF timeseries and netcdf depth-time grids,
using python and based on {ref}`xarray`.  The NetCDF files should be largely CF compliant.

The basic workflow is three glider-dependent steps to make a netcdf timeseries, followed
by optional steps to create profiles or depth-time grids.

```python
seaexplorer.raw_to_rawnc(rawdir, rawncdir, deploymentyaml)
seaexplorer.merge_rawnc(rawncdir, rawncdir, deploymentyaml, kind='sub')
outname = seaexplorer.raw_to_timeseries(rawncdir, tsdir, deploymentyaml, kind='sub')
# optional...
ncprocess.extract_timeseries_profiles(outname, profiledir, deploymentyaml)
outname2 = ncprocess.make_gridfiles(outname, griddir, deploymentyaml)
# optional plotting
pgplot.grid_plots(outname2, plottingyaml)
```

Data comes from and is written in directories, and metadata is supplied by a yaml file.

Currently only [Alseamar SeaExplorer](https://www.alseamar-alcen.com/products/underwater-glider/seaexplorer) and [Teledyne/Webb Slocum](http://www.teledynemarine.com/autonomous-underwater-gliders) glider data files are supported, and those with limited configurations.  Other gliders will hopefully be added.  If you have a glider type or configuration you would like added, [open an issue or pull request!](https://github.com/c-proof/pyglider).

```{toctree}
---
maxdepth: 1
---
Install
getting-started-seaexplorer
getting-started-slocum
plotting
pyglider/pyglider

```

## Acknowledgements

- Slocum binary translator based on
<https://gitlab.oceantrack.org/ocean-gliders-canada/dinkum/tree/seabBranch_py3>
- Processing steps closely follow the work by SOCIB
<https://github.com/socib/glider_toolbox>
- Rutger's description of the Slocum binary files is very helpful: <https://github.com/kerfoot/spt/wiki/Slocum-Glider-Data-File-Primer>
- The somewhat arcane metadata format for NGDAC is here: <https://github.com/ioos/ioosngdac/wiki/NGDAC-NetCDF-File-Format-Version-2>


Indices and tables
==================

* {ref}`genindex`
* {ref}`modindex`
* {ref}`search`
