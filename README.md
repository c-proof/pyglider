pyglider
========

Python tools for interacting with ocean glider data

Installation:
=============

- install conda/miniconda
- clone the repository onto your computer and cd into the directory
- from inside the cloned repository directory run  
`conda env create -f environment.yml`, it will read environment.yml and
automatically install the dependencies, and create environment `pyglider`,
(though you can change the name by editing `environment.yml`)

Editable installation
=====================

If you want to be able to edit the files in `pyglider/pyglider` then install
as above, and then do `pip install -e .`.  That will re-install pyglider
with links to the local directory, so you can edit the library files.  
If you do so consider making a pull-request with your changes!

Getting Started
===============

This repository has a fully-working example in `pyglider/example`.  Running
the script `process_deployment.py` should be able to process the raw data in
`example/binary` and make new versions of the files in the other
subdirectories.   

The work flow (for now) for a new deployment data set from a Slocum glider is

1. offload the data from the glider computer or get from the dockserver
2. copy data from `offload/Science/SENTLOGS/*.EBD`,
  `offload/Science/LOGS/*.EBD`, `offload/Main_board/SENTLOGS/*.DBD`,
  `offload/Science/LOGS/*.DBD`, to `./binary/` of your new deployment
  presumably in a different directory than `example`
3. copy the data from the caches to `yourdeployment/cac/`; these are
  found in `offload/Science/STATE/CACHE/*.CAC` and
  `offload/Main_board/STATE/CACHE/*.CAC`.
4. Edit `rosie_713_sensors.txt` (and probably name after your own glider)
  for the data streams from the raw files you want.  The name of all the
  variables are in the `*.CAC` files.  
5. Edit `deployment.yaml` with appropriate metadata for your glider and
  deployment.  There are a lot of fields, but most should not change after
  the first deploymnt unless the snsor payload changes on the glider.  

Why?
====

There are other glider workflows, but the only ones I know of are in Matlab,
and there is no reason for that.  

The only disadvantage of this code is that the dcoding step from binary to
the raw netcdf is quite slow for the Slocum binary data.  That can/should
be improved by using a C library instead of python, but really, its not *that*
slow.  
