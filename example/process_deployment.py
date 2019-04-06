import logging
import os
import pyglider.slocum as slocum
import pyglider.ncprocess as ncprocess
import pyglider.plotting as pgplot

logging.basicConfig(level='INFO')

binarydir  = './binary/'
rawdir     = './rawnc/'
cacdir     = './cac/'
sensorlist = './rosie_713_sensors.txt'
deploymentyaml = './deployment.yml'
l1tsdir    = './L1-timeseries/'
profiledir = './L1-profiles/'
griddir    = './L2-gridfiles/'
plottingyaml = './plottingconfig.yml'

if 0:
    # turn *.EBD and *.DBD into *.ebd.nc and *.dbd.nc netcdf files.
    slocum.binary_to_rawnc(binarydir, rawdir, cacdir, sensorlist, deploymentyaml)

    # merge individual neetcdf files into single netcdf files *.ebd.nc and *.dbd.nc
    slocum.merge_rawnc(rawdir, rawdir, deploymentyaml)

    # Make level-1 timeseries netcdf file from th raw files...
    outname = slocum.raw_to_L1timeseries(rawdir, l1tsdir, deploymentyaml)

    # make profile netcdf files for ioos gdac...
    ncprocess.extract_L1timeseries_profiles(outname, profiledir, deploymentyaml)

    # make grid of dataset....

    outname = 'L1-timeseries/rosie713-20180915T1202_L1.nc'
    pgplot.timeseries_plots(outname, plottingyaml)
outname = 'L1-timeseries/rosie713-20180915T1202_L1.nc'
outnamegrid = ncprocess.make_L2_gridfiles(outname, griddir, deploymentyaml)
outname = 'L2-gridfiles/rosie713-20180915T1202_L2grid.nc'
pgplot.grid_plots(outname, plottingyaml)
