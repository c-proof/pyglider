import logging
import os
import pyglider.seaexplorer as seaexplorer
import pyglider.ncprocess as ncprocess
import pyglider.plotting as pgplot

logging.basicConfig(level='INFO')

rawdir = './realtime_raw/'
rawncdir = './realtime_rawnc/'
deploymentyaml = './deploymentRealtime.yml'
l0tsdir = './L0-timeseries/'
profiledir = './L0-profiles/'
griddir = './L0-gridfiles/'
plottingyaml = './plottingconfig.yml'

if __name__ == '__main__':
    # clean last processing...
    os.system('rm ' + rawncdir + '* ' + l0tsdir + '* ' + profiledir + '* ' +
              griddir + '* ')

    # turn *.EBD and *.DBD into *.ebd.nc and *.dbd.nc netcdf files.
    seaexplorer.raw_to_rawnc(rawdir, rawncdir, deploymentyaml)
    # merge individual neetcdf files into single netcdf files *.ebd.nc and *.dbd.nc
    seaexplorer.merge_rawnc(rawncdir, rawncdir, deploymentyaml, kind='sub')
    # Make level-1 timeseries netcdf file from th raw files...
    outname = seaexplorer.raw_to_L0timeseries(rawncdir, l0tsdir, deploymentyaml, kind='sub')
    ncprocess.extract_L0timeseries_profiles(outname, profiledir, deploymentyaml)
    outname2 = ncprocess.make_L0_gridfiles(outname, griddir, deploymentyaml)
    # make profile netcdf files for ioos gdac...
    # make grid of dataset....
    pgplot.timeseries_plots(outname, plottingyaml)
    pgplot.grid_plots(outname2, plottingyaml)