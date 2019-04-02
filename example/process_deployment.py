import logging
import os
import pyglider.slocum as slocum
import pyglider.ncprocess as ncprocess

logging.basicConfig(level='DEBUG')

binarydir = './binary/'
rawdir    = './rawnc/'
cacdir    = './cac/'
sensorlist = './rosie_713_sensors.txt'
deploymentyaml = './deployment.yml'
l1tsdir    = './L1-timeseries/'
profiledir    = './L1-profiles/'


# turn *.EBD and *.DBD into *.ebd.nc and *.dbd.nc netcdf files.
slocum.binary_to_rawnc(binarydir, rawdir, cacdir, sensorlist, deploymentyaml)

# merge individual neetcdf files into single netcdf files *.ebd.nc and *.dbd.nc
slocum.merge_rawnc(rawdir, rawdir, deploymentyaml)

# Make level-1 timeseries netcdf file from th raw files...
outname = slocum.raw_to_L1timeseries(rawdir, l1tsdir, deploymentyaml)

# make profile netcdf files for ioos gdac...
ncprocess.extract_L1timeseries_profiles(outname, profiledir, deploymentyaml)
