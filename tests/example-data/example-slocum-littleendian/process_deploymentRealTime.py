import logging
import os
import shutil
import pyglider.ncprocess as ncprocess
import pyglider.utils as utils
import cproofutils.plotting as cpplot
import cproofutils.plotlinepmap as lpmap
import pyglider.slocum as slocum
import xarray as xr
import numpy as np
import locale

logging.basicConfig(level='INFO')

binarydir  = './realtime_raw/'
rawdir     = './realtime_rawnc/'
cacdir     = './cac/'
sensorlist = './dfo-maria997_sensors.txt'
deploymentyaml = './deployment.yml'
l1tsdir    = './L0-timeseries/'
profiledir = './L0-profiles/'
griddir    = './L0-gridfiles/'
plottingyaml = './plottingconfig.yml'
scisuffix    = 'tbd'
glidersuffix = 'sbd'

if False:
    os.system('rsync -av --no-perms --chmod=a+rX cproof@sfmc.webbresearch.com:/var/opt/sfmc-dockserver/stations/dfo/gliders/ ~/processing/slocum_dockserver/')
    # os.system('source synctodfo.sh')
    os.system('rsync -av ~/processing/slocum_dockserver/maria-997/from-glider/*2022-1* ' + binarydir)
    os.system('rsync -av ~/processing/slocum_dockserver/maria-997/from-glider/*2022-2* ' + binarydir)
    os.system('rsync -av ~/processing/slocum_dockserver/maria-997/gliderState.xml ./')
    for months in range(6, 6):
        os.system(f'rsync -av ~/processing/slocum_dockserver/maria-997/logs/maria-997_2022{months:02d}*.log ./logs/')


if True:
    # turn *.EBD and *.DBD into *.ebd.nc and *.dbd.nc netcdf files.
    slocum.binary_to_rawnc(binarydir, rawdir, cacdir, sensorlist, deploymentyaml,
            incremental=True, scisuffix=scisuffix, glidersuffix=glidersuffix)

    # remove some bad files:
    # os.system('rm realtime_rawnc/017*')
    # merge individual netcdf files into single netcdf files *.ebd.nc and *.dbd.nc
    slocum.merge_rawnc(rawdir, rawdir, deploymentyaml,
            scisuffix=scisuffix, glidersuffix=glidersuffix)
print('Done merge')

if True:
    # Make level-1 timeseries netcdf file from the raw files...
    outname = slocum.raw_to_timeseries(rawdir, l1tsdir, deploymentyaml)
    outname = utils.get_profiles(outname, filt_time=400, profile_min_time=100)
    # make profile netcdf files for ioos gdac...
    #ncprocess.extract_L1timeseries_profiles(outname, profiledir, deploymentyaml)
    # make grid of dataset....


if False:
    cpplot.timeseries_plots(outname, plottingyaml)
    # correct the timeseries file if we don't have sbd data for lat,lon, time

if True:
    outname2 = ncprocess.make_gridfiles(outname, griddir, deploymentyaml,
                                           dz=10)
    cpplot.grid_plots(outname2, plottingyaml)
