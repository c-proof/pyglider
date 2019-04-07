import glob
import logging
import os
import pyglider.slocum as slocum
import pyglider.ncprocess as ncprocess

logging.basicConfig(level='DEBUG')

binarydir  = './binary/'
cacdir  = './cac/'

# turn *.EBD and *.DBD into *.ebd.nc and *.dbd.nc netcdf files.
d = binarydir + '*.EBD'
filesScience = glob.glob(d)
filesScience.sort()

for fn in filesScience:
    try:
        fmeta, _ = slocum.dbd_get_meta(fn, cachedir=cacdir)
        print(fn, fmeta['fileopen_time'])
    except:
        print('Could not read ', fn)
        pass
