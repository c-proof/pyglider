# -*- coding: utf-8 -*-
"""
Process dfo-rosie713 realtime data to OG 1.0 format.

This script mirrors tests/example-slocum/process_deploymentRealTime.py but
targets OG 1.0 output.  It uses the same raw binary data; only the YAML and
output directories differ.

NOTE: This script documents the *intended* API once pyglider has been updated
to support processing_role, processing_method, output_dimension, and the
_load_dataset/_save_dataset helpers.  Some calls will not work correctly until
those changes are made to slocum.py, utils.py, and ncprocess.py.
"""

import logging
import os
import pyglider.ncprocess as ncprocess
import pyglider.slocum as slocum
import pyglider.utils as pgutils

logging.basicConfig(level='INFO')

# Raw data lives in the existing example-slocum directory.
# Outputs go into subdirectories here.
binarydir = '../example-slocum/realtime_raw/'
cacdir    = '../example-slocum/cac/'

deploymentyaml = './deploymentRealtime_og10.yml'

l1tsdir  = './L0-timeseries/'
profiledir = './L0-profiles/'
griddir  = './L0-gridfiles/'

# ------------------------------------------------------------------------
# Step 1: binary → OG 1.0 timeseries netCDF
#
# binary_to_timeseries will need to:
#   - read processing_role to identify pressure/temperature/conductivity/
#     latitude/longitude variables by role rather than by hardcoded name
#   - read processing_method entries to compute derived variables
#     (DEPTH, PSAL, SIGMA0, DENSITY, POTENTIAL_TEMPERATURE, PROFILE_NUMBER,
#      PROFILE_DIRECTION, DISTANCE_OVER_GROUND) using the named inputs
#   - call _save_dataset instead of ds.to_netcdf so that the time dimension
#     is renamed to N_MEASUREMENTS before writing
# ------------------------------------------------------------------------
outname = slocum.binary_to_timeseries(
    binarydir, cacdir, l1tsdir, deploymentyaml,
    search='*.[s|t]bd',
    profile_filt_time=20,
    profile_min_time=20,
)

# ------------------------------------------------------------------------
# Step 2: timeseries → per-profile netCDF files (IOOS GDAC style)
#
# extract_timeseries_profiles will need to:
#   - call _load_dataset instead of xr.open_dataset so that the
#     N_MEASUREMENTS dimension is normalised back to time for processing
#   - use processing_role to find profile_index (PROFILE_NUMBER),
#     latitude (LATITUDE), longitude (LONGITUDE), etc.
#   - call _save_dataset when writing each profile file
#
# NOTE: OG 1.0 does not define per-profile files; the trajectory file is
# the primary deliverable.  This step produces IOOS GDAC profile files as
# a secondary output for users who need them.
# ------------------------------------------------------------------------
# ncprocess.extract_timeseries_profiles(outname, profiledir, deploymentyaml)

# ------------------------------------------------------------------------
# Step 3: timeseries → gridded netCDF
#
# make_gridfiles will need to:
#   - call _load_dataset to normalise the dimension on load
#   - use processing_role to find depth (DEPTH), latitude (LATITUDE),
#     longitude (LONGITUDE) for gridding axes
#   - call _save_dataset when writing the gridded file
# ------------------------------------------------------------------------
outname2 = ncprocess.make_gridfiles(outname, griddir, deploymentyaml)

pgutils.example_gridplot(
    outname2,
    './gridplot_og10.png',
    ylim=[150, 0],
    toplot=['POTENTIAL_TEMPERATURE', 'PSAL', 'DOXY', 'CHLA'],
)
