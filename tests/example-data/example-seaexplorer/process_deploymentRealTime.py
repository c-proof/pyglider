import logging
import os
import pyglider.seaexplorer as seaexplorer
import pyglider.ncprocess as ncprocess
import pyglider.plotting as pgplot

logging.basicConfig(level='INFO')

sourcedir = '~alseamar/Documents/SEA035/000012/000012/C-Csv/*'
rawdir  = './realtime_raw/'
rawncdir     = './realtime_rawnc/'
deploymentyaml = './deploymentRealtime.yml'
l0tsdir    = './L0-timeseries/'
profiledir = './L0-profiles/'
griddir    = './L0-gridfiles/'
plottingyaml = './plottingconfig.yml'

## get the data and clean up derived
#os.system('source synctodfo.sh')
if 0:
    os.system('rsync -av ' + sourcedir + ' ' + rawdir)

# clean last processing...
os.system('rm ' + rawncdir + '* ' + l0tsdir + '* ' + profiledir + '* ' +
          griddir + '* ')


if 1:
    # turn *.EBD and *.DBD into *.ebd.nc and *.dbd.nc netcdf files.
    seaexplorer.raw_to_rawnc(rawdir, rawncdir, deploymentyaml)
        # merge individual neetcdf files into single netcdf files *.ebd.nc and *.dbd.nc
    seaexplorer.merge_rawnc(rawncdir, rawncdir, deploymentyaml, kind='sub')

        # Make level-1 timeseries netcdf file from th raw files...
    outname = seaexplorer.raw_to_timeseries(rawncdir, l0tsdir, deploymentyaml, kind='sub')
    ncprocess.extract_timeseries_profiles(outname, profiledir, deploymentyaml)
    outname2 = ncprocess.make_gridfiles(outname, griddir, deploymentyaml)

if 1:
    # make profile netcdf files for ioos gdac...
    # make grid of dataset....
    pgplot.timeseries_plots(outname, plottingyaml)
    pgplot.grid_plots(outname2, plottingyaml)
