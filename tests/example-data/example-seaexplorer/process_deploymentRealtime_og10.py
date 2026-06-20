# -*- coding: utf-8 -*-
"""
Process dfo-eva035 realtime data to OG 1.0 format.

This script mirrors process_deploymentRealTime.py but targets OG 1.0 output.
It reuses the raw parquet files produced by the standard pipeline; only the
YAML and output directories differ.
"""

import logging
import pyglider.seaexplorer as seaexplorer
import pyglider.ncprocess as ncprocess

logging.basicConfig(level='INFO')

rawncdir       = './realtime_rawnc/'
deploymentyaml = './deploymentRealtime_og10.yml'
l0tsdir        = './L0-timeseries-og10/'
griddir        = './L0-gridfiles-og10/'

# Step 1: timeseries (raw_to_rawnc and merge_parquet already run by the
# standard pipeline — we reuse the parquet files in realtime_rawnc/).
outname = seaexplorer.raw_to_timeseries(
    rawncdir, l0tsdir, deploymentyaml, kind='sub'
)

# Step 2: gridded netCDF
outname2 = ncprocess.make_gridfiles(outname, griddir, deploymentyaml)
