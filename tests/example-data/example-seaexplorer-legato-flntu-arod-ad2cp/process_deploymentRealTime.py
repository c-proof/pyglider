import logging
import os
import pyglider.seaexplorer as seaexplorer
import pyglider.ncprocess as ncprocess
import pyglider.utils as pgutils


logging.basicConfig(level='INFO')

rawdir = './realtime_raw/'
rawncdir = './realtime_rawnc/'
deploymentyaml = './deploymentRealtime.yml'
l0tsdir = './L0-timeseries/'
profiledir = './L0-profiles/'
griddir = './L0-gridfiles/'

if __name__ == '__main__':
    # clean last processing...
    os.system('rm ' + rawncdir + '* ' + l0tsdir + '* ' + profiledir + '* ' +
              griddir + '* ')

    # turn seaexplorer zipped csvs into nc files.
    seaexplorer.raw_to_rawnc(rawdir, rawncdir, deploymentyaml)
    # merge individual netcdf files into single netcdf files *.gli*.nc and *.pld1*.nc
    seaexplorer.merge_parquet(rawncdir, rawncdir, deploymentyaml, kind='sub')
    # Make level-0 timeseries netcdf file from the raw files...
    outname = seaexplorer.raw_to_timeseries(rawncdir, l0tsdir, deploymentyaml, kind='sub')
    outname = pgutils.get_profiles(outname)
    ncprocess.extract_timeseries_profiles(outname, profiledir, deploymentyaml)
    outname2 = ncprocess.make_gridfiles(outname, griddir, deploymentyaml)
    # make profile netcdf files for ioos gdac...
    # make grid of dataset....
    pgutils.example_gridplot(outname2, './gridplot.png', ylim=[100, 0],
                             toplot=['potential_temperature', 'salinity', 'oxygen_concentration',
                                     'chlorophyll', 'cdom'])
