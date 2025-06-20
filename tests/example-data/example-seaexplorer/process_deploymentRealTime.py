import logging
import os
import pyglider.seaexplorer as seaexplorer
import pyglider.ncprocess as ncprocess
import pyglider.utils as pgutils

logging.basicConfig(level='INFO')

sourcedir = '~alseamar/Documents/SEA035/000012/000012/C-Csv/*'
rawdir  = './realtime_raw/'
rawncdir     = './realtime_rawnc/'
deploymentyaml = './deploymentRealtime.yml'
l0tsdir    = './L0-timeseries/'
profiledir = './L0-profiles/'
griddir    = './L0-gridfiles/'

## get the data and clean up derived
if False:
    os.system('rsync -av ' + sourcedir + ' ' + rawdir)

# clean last processing...
os.system('rm ' + rawncdir + '* ' + l0tsdir + '* ' + profiledir + '* ' +
          griddir + '* ')

if True:
    # turn *.EBD and *.DBD into *.ebd.nc and *.dbd.nc netcdf files.
    seaexplorer.raw_to_rawnc(rawdir, rawncdir, deploymentyaml)
        # merge individual neetcdf files into single netcdf files *.ebd.nc and *.dbd.nc
    seaexplorer.merge_parquet(rawncdir, rawncdir, deploymentyaml, kind='sub')

        # Make level-1 timeseries netcdf file from th raw files...
    outname = seaexplorer.raw_to_timeseries(rawncdir, l0tsdir, deploymentyaml, kind='sub',
                                            deadreckon=False)
    ncprocess.extract_timeseries_profiles(outname, profiledir, deploymentyaml)
    outname2 = ncprocess.make_gridfiles(outname, griddir, deploymentyaml)

    pgutils.example_gridplot(outname2, './gridplot.png', ylim=[700, 0],
                             toplot=['potential_temperature', 'salinity', 'oxygen_concentration',
                                     'chlorophyll', 'cdom'])
