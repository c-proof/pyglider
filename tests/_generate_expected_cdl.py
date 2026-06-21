"""
Regenerate NC golden files and CDL headers for the CDL-based regression tests.

Run this after an intentional change to update the expected values:

    cd <repo root>
    python tests/_generate_expected_cdl.py

Each golden file pair consists of:
  - A cleaned .nc file (SKIP_ATTRS stripped) stored in tests/expected/
  - A .cdl header file (ncdump -h output) alongside it for git diffs

Tests compare pipeline output against the golden .nc files using
xarray.testing.assert_identical (via nc_test_helpers.assert_datasets_equal),
after stripping dynamic attributes from the actual output.

The script re-runs the seaexplorer and process_adjusted pipelines itself so
that the golden files always reflect the current code.  The slocum L0
timeseries pipeline requires dbdreader (a binary reader not available here),
so the slocum golden file is generated from the NC file produced by
process_deploymentRealTime.py in tests/example-data/example-slocum/; run
that script before running this one if the slocum golden file needs updating.
"""

import shutil
import subprocess
from pathlib import Path

import xarray as xr

import pyglider.seaexplorer as seaexplorer
import pyglider.ncprocess as ncprocess
import pyglider.slocum as slocum
from pyglider.process_adjusted import run_process_adjusted


SKIP_ATTRS = {'date_created', 'date_issued', 'history'}


def save_golden(nc_path: Path, golden_nc_path: Path) -> None:
    """Strip dynamic attrs, write cleaned NC golden file, and generate CDL header.

    If *nc_path* and *golden_nc_path* are the same file, the file is rewritten
    in place via a temporary path.
    """
    nc_path = Path(nc_path)
    golden_nc_path = Path(golden_nc_path)
    golden_nc_path.parent.mkdir(parents=True, exist_ok=True)

    with xr.open_dataset(nc_path) as ds:
        ds = ds.load()  # pull into memory before closing
        for attr in SKIP_ATTRS:
            ds.attrs.pop(attr, None)
        if golden_nc_path.resolve() == nc_path.resolve():
            tmp = nc_path.with_suffix('.tmp.nc')
            ds.to_netcdf(tmp)
            tmp.replace(nc_path)
        else:
            ds.to_netcdf(golden_nc_path)

    cdl_path = golden_nc_path.with_suffix('.cdl')
    with open(cdl_path, 'w') as f:
        subprocess.run(['ncdump', '-h', str(golden_nc_path)], stdout=f, check=True)
    print(f'wrote {golden_nc_path} and {cdl_path}')


def remove_old_yaml(expected_dir: Path) -> None:
    """Remove stale YAML golden files left over from the previous YAML-based system."""
    removed = list(expected_dir.rglob('*.yml'))
    for p in removed:
        p.unlink()
        print(f'removed {p}')


if __name__ == '__main__':
    tests_dir = Path(__file__).parent
    expected_dir = tests_dir / 'expected'
    example_data = tests_dir / 'example-data'

    # ------------------------------------------------------------------
    # Slocum OG 1.0 pipeline
    # ------------------------------------------------------------------
    sl_og10_cacdir = example_data / 'example-slocum/cac/'
    sl_og10_binary = str(example_data / 'example-slocum/realtime_raw/') + '/'
    sl_og10_yaml = str(example_data / 'example-slocum/deploymentRealtime_og10.yml')
    sl_og10_ts = str(example_data / 'example-slocum/L0-timeseries-og10/') + '/'
    sl_og10_nc = slocum.binary_to_timeseries(
        sl_og10_binary, sl_og10_cacdir, sl_og10_ts, sl_og10_yaml,
        search='*.[s|t]bd', profile_filt_time=20, profile_min_time=20,
    )

    # ------------------------------------------------------------------
    # SeaExplorer NRT sub pipeline
    # ------------------------------------------------------------------
    se_raw = str(example_data / 'example-seaexplorer/realtime_raw/') + '/'
    se_rawnc = str(example_data / 'example-seaexplorer/realtime_rawnc/') + '/'
    se_yaml = str(example_data / 'example-seaexplorer/deploymentRealtime.yml')
    se_l0ts = str(example_data / 'example-seaexplorer/L0-timeseries-test/') + '/'
    se_l0ts_adj = str(example_data / 'example-seaexplorer/L0-timeseries/') + '/'
    seaexplorer.raw_to_rawnc(se_raw, se_rawnc, se_yaml)
    seaexplorer.merge_parquet(se_rawnc, se_rawnc, se_yaml, kind='sub')
    se_nc = seaexplorer.raw_to_L0timeseries(se_rawnc, se_l0ts, se_yaml, kind='sub')
    # Also write to L0-timeseries/ so run_process_adjusted has fresh input
    seaexplorer.raw_to_L0timeseries(se_rawnc, se_l0ts_adj, se_yaml, kind='sub')

    # ------------------------------------------------------------------
    # SeaExplorer OG 1.0 pipeline (reuses parquet files from NRT sub above)
    # ------------------------------------------------------------------
    se_og10_yaml = str(example_data / 'example-seaexplorer/deploymentRealtime_og10.yml')
    se_og10_ts = str(example_data / 'example-seaexplorer/L0-timeseries-og10/') + '/'
    se_og10_griddir = str(example_data / 'example-seaexplorer/L0-gridfiles-og10/') + '/'
    se_og10_nc = seaexplorer.raw_to_timeseries(se_rawnc, se_og10_ts, se_og10_yaml, kind='sub')
    se_og10_grid_nc = ncprocess.make_gridfiles(se_og10_nc, se_og10_griddir, se_og10_yaml)

    # ------------------------------------------------------------------
    # SeaExplorer raw delayed pipeline
    # ------------------------------------------------------------------
    ser_raw = str(example_data / 'example-seaexplorer-raw/delayed_raw/') + '/'
    ser_rawnc = str(example_data / 'example-seaexplorer-raw/delayed_rawnc/') + '/'
    ser_yaml = str(example_data / 'example-seaexplorer-raw/deployment.yml')
    ser_l0ts = str(example_data / 'example-seaexplorer-raw/L0-timeseries-test/') + '/'
    seaexplorer.raw_to_rawnc(ser_raw, ser_rawnc, ser_yaml)
    seaexplorer.merge_parquet(ser_rawnc, ser_rawnc, ser_yaml, kind='raw')
    ser_nc = seaexplorer.raw_to_L0timeseries(ser_rawnc, ser_l0ts, ser_yaml, kind='raw')

    # ------------------------------------------------------------------
    # process_adjusted pipelines — write directly to expected_dir so the
    # output NC files become the golden files (cleaned in place below)
    # ------------------------------------------------------------------
    se_adj_nc = run_process_adjusted(
        expected_dir / 'example-seaexplorer',
        deploy_name='dfo-eva035-20190718',
        deployfile=example_data / 'example-seaexplorer/deploymentRealtime.yml',
        adjustedyaml=example_data / 'example-seaexplorer/adjusted.yml',
        input_dir=example_data / 'example-seaexplorer',
    )
    sl_adj_nc = run_process_adjusted(
        expected_dir / 'example-slocum',
        deploy_name='dfo-rosie713-20190615',
        deployfile=example_data / 'example-slocum/deploymentRealtime.yml',
        adjustedyaml=example_data / 'example-slocum/adjusted.yml',
        input_dir=example_data / 'example-slocum',
    )

    # ------------------------------------------------------------------
    # Map NC source files to golden NC destinations.
    # For process_adjusted outputs the src == dst (cleaned in place).
    # ------------------------------------------------------------------
    todo = {
        # slocum L0 timeseries — requires dbdreader; read pre-existing NC
        example_data / 'example-slocum/L0-timeseries/dfo-rosie713-20190615.nc':
            expected_dir / 'example-slocum/L0-timeseries/dfo-rosie713-20190615.nc',

        # slocum OG 1.0
        Path(sl_og10_nc):
            expected_dir / 'example-slocum/L0-timeseries-og10/dfo-rosie713-20190615.nc',

        # seaexplorer NRT sub
        Path(se_nc):
            expected_dir / 'example-seaexplorer/L0-timeseries/dfo-eva035-20190718.nc',

        # seaexplorer OG 1.0 timeseries
        Path(se_og10_nc):
            expected_dir / 'example-seaexplorer/L0-timeseries-og10/dfo-eva035-20190718.nc',

        # seaexplorer OG 1.0 gridfile
        Path(se_og10_grid_nc):
            expected_dir / 'example-seaexplorer/L0-gridfiles-og10/dfo-eva035-20190718.nc',

        # seaexplorer raw delayed
        Path(ser_nc):
            expected_dir / 'example-seaexplorer-raw/L0-timeseries/dfo-bb046-20200908.nc',

        # process_adjusted outputs (src == dst, cleaned in place)
        Path(se_adj_nc):
            Path(se_adj_nc),
        expected_dir / 'example-seaexplorer/L0-gridfiles/dfo-eva035-20190718_grid_adjusted.nc':
            expected_dir / 'example-seaexplorer/L0-gridfiles/dfo-eva035-20190718_grid_adjusted.nc',
        Path(sl_adj_nc):
            Path(sl_adj_nc),
        expected_dir / 'example-slocum/L0-gridfiles/dfo-rosie713-20190615_grid_adjusted.nc':
            expected_dir / 'example-slocum/L0-gridfiles/dfo-rosie713-20190615_grid_adjusted.nc',
    }

    for nc_src, nc_dst in todo.items():
        if not nc_src.exists():
            print(f'SKIP (NC not found): {nc_src}')
            continue
        save_golden(nc_src, nc_dst)

    # Remove stale YAML golden files from the old system
    remove_old_yaml(expected_dir)
