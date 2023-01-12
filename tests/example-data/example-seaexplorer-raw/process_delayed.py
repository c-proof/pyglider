import logging
import os
import pyglider.seaexplorer as seaexplorer

logging.basicConfig(level='INFO')

rawdir  = './delayed_raw/'
rawncdir     = './realtime_rawnc/'
deploymentyaml = './deployment.yml'
l0tsdir    = './L0-timeseries/'
profiledir = './L0-profiles/'
griddir    = './L0-gridfiles/'


if __name__ == '__main__':
    # clean last processing
    os.system('rm ' + rawncdir + '* ' + l0tsdir + '* ' + profiledir + '* ' +
              griddir + '* ')
    # Convert from csv to parquet
    seaexplorer.raw_to_rawnc(rawdir, rawncdir, deploymentyaml)
    # merge individual parquet files
    seaexplorer.merge_parquet(rawncdir, rawncdir, deploymentyaml, kind='raw')

    # Make level-1 timeseries netcdf file from the raw files
    outname = seaexplorer.raw_to_timeseries(rawncdir, l0tsdir, deploymentyaml, kind='raw')
