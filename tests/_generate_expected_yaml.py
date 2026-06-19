"""
Regenerate YAML summary files for the YAML-based regression tests.

Run this after an intentional change to update the expected values:

    cd <repo root>
    python tests/_generate_expected_yaml.py

Each YAML file stores per-variable statistics (min, max, mean, std, n_valid,
n_nan) and time-derivative statistics (diff_mean, diff_std), plus global and
variable-level CF attributes.  These are compared by test_*_yaml.py tests.

The script re-runs the seaexplorer and process_adjusted pipelines itself so
that the golden files always reflect the current code.  The slocum L0
timeseries pipeline requires dbdreader (a binary reader not available here),
so that golden file is generated from the NC file produced by
process_deploymentRealTime.py in tests/example-data/example-slocum/; run
that script before running this one if the slocum golden file needs updating.
"""

from pathlib import Path

import numpy as np
import xarray as xr
import yaml

import pyglider.seaexplorer as seaexplorer
from pyglider.process_adjusted import run_process_adjusted


SKIP_ATTRS = {'date_created', 'date_issued', 'history'}


def summarize(ds: xr.Dataset) -> dict:
    result = {
        'attrs': {k: str(v) for k, v in ds.attrs.items()
                  if k not in SKIP_ATTRS},
        'variables': {},
    }
    for var in ds.variables:
        da = ds[var]
        vals = da.values

        if np.issubdtype(vals.dtype, np.datetime64):
            vals = vals.astype('datetime64[ns]').astype('float64')

        vals = vals.flatten().astype('float64')
        valid = vals[~np.isnan(vals)]
        diffs = np.diff(valid)

        entry = {
            'attrs': {k: str(v) for k, v in da.attrs.items()},
            'n_valid': int(len(valid)),
            'n_nan': int(np.sum(np.isnan(vals))),
        }
        if len(valid):
            entry.update(
                min=float(f'{np.min(valid):.12g}'),
                max=float(f'{np.max(valid):.12g}'),
                mean=float(f'{np.mean(valid):.12g}'),
                std=float(f'{np.std(valid):.12g}'),
            )
        if len(diffs):
            entry.update(
                diff_mean=float(f'{np.mean(diffs):.12g}'),
                diff_std=float(f'{np.std(diffs):.12g}'),
            )
        result['variables'][var] = entry
    return result


def write(nc_path: Path, yaml_path: Path) -> None:
    yaml_path.parent.mkdir(parents=True, exist_ok=True)
    ds = xr.open_dataset(nc_path)
    summary = summarize(ds)
    with open(yaml_path, 'w') as f:
        yaml.dump(summary, f, default_flow_style=False, allow_unicode=True,
                  sort_keys=True)
    print(f'wrote {yaml_path}')


if __name__ == '__main__':
    tests_dir = Path(__file__).parent
    expected_yaml = tests_dir / 'expected'
    example_data = tests_dir / 'example-data'

    # ------------------------------------------------------------------
    # SeaExplorer NRT sub pipeline
    # ------------------------------------------------------------------
    se_raw = str(example_data / 'example-seaexplorer/realtime_raw/') + '/'
    se_rawnc = str(example_data / 'example-seaexplorer/realtime_rawnc/') + '/'
    se_yaml = str(example_data / 'example-seaexplorer/deploymentRealtime.yml')
    se_l0ts = str(example_data / 'example-seaexplorer/L0-timeseries-test/') + '/'
    seaexplorer.raw_to_rawnc(se_raw, se_rawnc, se_yaml)
    seaexplorer.merge_parquet(se_rawnc, se_rawnc, se_yaml, kind='sub')
    se_nc = seaexplorer.raw_to_L0timeseries(se_rawnc, se_l0ts, se_yaml, kind='sub')

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
    # process_adjusted pipelines (seaexplorer and slocum)
    # ------------------------------------------------------------------
    se_adj_nc = run_process_adjusted(
        expected_yaml / 'example-seaexplorer',
        deploy_name='dfo-eva035-20190718',
        deployfile=example_data / 'example-seaexplorer/deploymentRealtime.yml',
        adjustedyaml=example_data / 'example-seaexplorer/adjusted.yml',
        input_dir=example_data / 'example-seaexplorer',
    )
    sl_adj_nc = run_process_adjusted(
        expected_yaml / 'example-slocum',
        deploy_name='dfo-rosie713-20190615',
        deployfile=example_data / 'example-slocum/deploymentRealtime.yml',
        adjustedyaml=example_data / 'example-slocum/adjusted.yml',
        input_dir=example_data / 'example-slocum',
    )

    # ------------------------------------------------------------------
    # Map NC files to YAML golden files
    # ------------------------------------------------------------------
    todo = {
        # slocum L0 timeseries — requires dbdreader, read from existing NC
        example_data / 'example-slocum/L0-timeseries/dfo-rosie713-20190615.nc':
            expected_yaml / 'example-slocum/L0-timeseries/dfo-rosie713-20190615.yml',

        # seaexplorer NRT sub — just generated above
        Path(se_nc):
            expected_yaml / 'example-seaexplorer/L0-timeseries/dfo-eva035-20190718.yml',

        # seaexplorer raw delayed — just generated above
        Path(ser_nc):
            expected_yaml / 'example-seaexplorer-raw/L0-timeseries/dfo-bb046-20200908.yml',

        # process_adjusted outputs — just generated above
        Path(se_adj_nc):
            expected_yaml / 'example-seaexplorer/L0-timeseries/dfo-eva035-20190718_adjusted.yml',
        expected_yaml / 'example-seaexplorer/L0-gridfiles/dfo-eva035-20190718_grid_adjusted.nc':
            expected_yaml / 'example-seaexplorer/L0-gridfiles/dfo-eva035-20190718_grid_adjusted.yml',
        Path(sl_adj_nc):
            expected_yaml / 'example-slocum/L0-timeseries/dfo-rosie713-20190615_adjusted.yml',
        expected_yaml / 'example-slocum/L0-gridfiles/dfo-rosie713-20190615_grid_adjusted.nc':
            expected_yaml / 'example-slocum/L0-gridfiles/dfo-rosie713-20190615_grid_adjusted.yml',
    }

    for nc, yml in todo.items():
        if not nc.exists():
            print(f'SKIP (NC not found): {nc}')
            continue
        write(nc, yml)
