# -*- coding: utf-8 -*-
import logging
import os
import pyglider.ncprocess as ncprocess
import pyglider.slocum as slocum
import pyglider.utils as pgutils

logging.basicConfig(level='INFO')

binarydir = './realtime_raw/'
rawdir = './realtime_rawnc/'
cacdir = './cac/'
sensorlist = './dfo-rosie713_sensors.txt'
deploymentyaml = './deploymentRealtime.yml'
l1tsdir = './L0-timeseries/'
profiledir = './L0-profiles/'
griddir = './L0-gridfiles/'
scisuffix = 'tbd'
glidersuffix = 'sbd'

do_direct = True
# only do this for a real run, or something like this
real = False
if real:
    os.system('rsync -av cproof@sfmc.webbresearch.com:/var/opt/sfmc-dockserver/' +
              'stations/dfo/gliders/ ~/gliderdata/slocum_dockserver/')
    os.system('rsync -av ~/gliderdata/slocum_dockserver/rosie_713/from-glider/* ' +
              binarydir)

    os.system('rm ' + rawdir + 'dfo* ' + rawdir + 'TEMP*.nc ' + l1tsdir + '* ' +
              profiledir + '* ' + griddir + '* ')

if do_direct:
    # turn *.sdb and *.tbd into timeseries netcdf files
    outname = slocum.binary_to_timeseries(
        binarydir, cacdir, l1tsdir, deploymentyaml, search='*.[s|t]bd',
        profile_filt_time=20, profile_min_time=20)
else:
    # turn *.EBD and *.DBD into *.ebd.nc and *.dbd.nc netcdf files.
    slocum.binary_to_rawnc(
        binarydir, rawdir, cacdir, sensorlist, deploymentyaml,
        incremental=True, scisuffix=scisuffix, glidersuffix=glidersuffix)

    # merge individual neetcdf files into single netcdf files *.ebd.nc and *.dbd.nc
    slocum.merge_rawnc(
        rawdir, rawdir, deploymentyaml,
        scisuffix=scisuffix, glidersuffix=glidersuffix)

    # Make level-1 timeseries netcdf file from th raw files...
    outname = slocum.raw_to_timeseries(
        rawdir, l1tsdir, deploymentyaml,
        profile_filt_time=100, profile_min_time=300)

if True:
    # make profile netcdf files for ioos gdac...
    ncprocess.extract_timeseries_profiles(outname, profiledir, deploymentyaml)

# make grid of dataset....

outname2 = ncprocess.make_gridfiles(outname, griddir, deploymentyaml)
pgutils.example_gridplot(outname2, './gridplot2.png', ylim=[150, 0],
                         toplot=['potential_temperature', 'salinity',
                                 'oxygen_concentration', 'chlorophyll', 'cdom'])
