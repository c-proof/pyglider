import logging
from pathlib import Path
import xarray as xr

import pyglider.ncprocess as ncprocess
import pyglider.utils as utils

logging.basicConfig(level=logging.INFO)
_log = logging.getLogger(__name__)


def run_process_adjusted(
    base_dir,
    deploy_name=None,
    deployfile=None,
    adjustedyaml=None,
):
    base_dir = Path(base_dir)

    deploy_name = deploy_name or base_dir.name

    # --- Paths ---
    ts_path = base_dir / 'L0-timeseries'
    gridpath = base_dir / 'L0-gridfiles'

    openfile = ts_path / f'{deploy_name}.nc'
    deployfile = Path(deployfile) if deployfile else base_dir / 'deploymentRealtime.yml'
    adjustedyaml = Path(adjustedyaml) if adjustedyaml else base_dir / 'adjusted.yml'

    # --- Load ---
    ts = xr.open_dataset(openfile)
    deployment = utils._get_deployment([deployfile, adjustedyaml])

    # --- Processing ---
    ts = utils.flag_CTD_data(ts)
    ts = utils.adjust_CTD(ts, deployment)

    # --- Save ---
    outfile = ts_path / f'{deploy_name}_adjusted.nc'
    _log.info('Saving adjusted timeseries to netcdf')
    ts.to_netcdf(outfile)

    # --- Grid ---
    ncprocess.make_gridfiles(
        str(outfile),
        str(gridpath),
        [str(deployfile), str(adjustedyaml)],
        fnamesuffix='_adjusted',
        maskfunction=utils.maskQC4,
    )

    return outfile


if __name__ == "__main__":
    run_process_adjusted(Path.cwd())