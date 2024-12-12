import logging
import os
import pyglider.seaexplorer as seaexplorer
logging.basicConfig(level='INFO')

rawdir = './realtime_raw/'
rawncdir = './realtime_rawnc/'
deploymentyaml = './ocean_gliders.yml'
l0tsdir = './L0-timeseries/'

if __name__ == '__main__':
    # clean last processing...
    os.system('rm ' + rawncdir + '* ' + l0tsdir + '* ')

    # turn seaexplorer zipped csvs into nc files.
    seaexplorer.raw_to_rawnc(rawdir, rawncdir, deploymentyaml)
    # merge individual netcdf files into single netcdf files *.gli*.nc and *.pld1*.nc
    seaexplorer.merge_parquet(rawncdir, rawncdir, deploymentyaml, kind='sub')
    # Make OceanGliders valid netcdf file
    outname = seaexplorer.raw_to_timeseries(rawncdir, l0tsdir, deploymentyaml, kind='sub', og_format=True)


