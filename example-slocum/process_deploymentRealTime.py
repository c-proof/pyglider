import logging
import os
import pyglider.ncprocess as ncprocess
import pyglider.plotting as pgplot
import pyglider.slocum as slocum

logging.basicConfig(level='INFO')

binarydir  = './realtime_raw/'
rawdir     = './realtime_rawnc/'
cacdir     = './cac/'
sensorlist = './dfo-rosie713_sensors.txt'
deploymentyaml = './deploymentRealtime.yml'
l1tsdir    = './L0-timeseries/'
profiledir = './L0-profiles/'
griddir    = './L0-gridfiles/'
plottingyaml = './plottingconfig.yml'
scisuffix    = 'tbd'
glidersuffix = 'sbd'


# only do this for a real run, or something like this
real = False
if real:
    os.system('rsync -av cproof@sfmc.webbresearch.com:/var/opt/sfmc-dockserver/stations/dfo/gliders/ ~/gliderdata/slocum_dockserver/')
    os.system('rsync -av ~/gliderdata/slocum_dockserver/rosie_713/from-glider/* ' + binarydir)

os.system('rm ' + rawdir + 'dfo* ' + rawdir + 'TEMP*.nc ' + l1tsdir + '* ' + profiledir + '* ' +
          griddir + '* ')

print(scisuffix)
# turn *.EBD and *.DBD into *.ebd.nc and *.dbd.nc netcdf files.
if 1:
    slocum.binary_to_rawnc(binarydir, rawdir, cacdir, sensorlist, deploymentyaml,
        incremental=True, scisuffix=scisuffix, glidersuffix=glidersuffix)

# merge individual neetcdf files into single netcdf files *.ebd.nc and *.dbd.nc
    slocum.merge_rawnc(rawdir, rawdir, deploymentyaml,
        scisuffix=scisuffix, glidersuffix=glidersuffix)
if 1:
    # Make level-1 timeseries netcdf file from th raw files...
    outname = slocum.raw_to_timeseries(rawdir, l1tsdir, deploymentyaml,
              profile_filt_time=100, profile_min_time=300)
    # make profile netcdf files for ioos gdac...
    ncprocess.extract_timeseries_profiles(outname, profiledir, deploymentyaml)

    # make grid of dataset....

outname = ncprocess.make_gridfiles(outname, griddir, deploymentyaml)
pgplot.grid_plots(outname, plottingyaml)
