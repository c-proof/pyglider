"""
YAML-summary-based regression tests for the example-slocum L0 timeseries.

Intended to replace test_slocum.py.  Instead of comparing full binary NetCDF
golden files, these tests check:
  - variable list
  - global attributes (excluding dynamic/version fields)
  - per-variable CF attributes
  - per-variable statistics (min, max, mean, std, n_valid, n_nan)
  - per-variable time-derivative statistics (diff_mean, diff_std)
  - time monotonicity
  - CF and GliderDAC compliance of profiles and timeseries

To regenerate the YAML after an intentional change, run:
    python tests/_generate_expected_yaml.py
"""

import json
from pathlib import Path

import numpy as np
import pytest
import xarray as xr
import yaml
from compliance_checker.runner import CheckSuite, ComplianceChecker

import pyglider.ncprocess as ncprocess
import pyglider.slocum as slocum

from yaml_test_helpers import RTOL, ATOL, SKIP_ATTRS, to_float, stats

library_dir = Path(__file__).parent.parent.absolute()
example_dir = library_dir / 'tests/example-data/'

# ---------------------------------------------------------------------------
# Run the pipeline once at module level
# ---------------------------------------------------------------------------
cacdir = example_dir / 'example-slocum/cac/'
binarydir = str(example_dir / 'example-slocum/realtime_raw/') + '/'
deploymentyaml = str(example_dir / 'example-slocum/deploymentRealtime.yml')
tsdir = str(example_dir / 'example-slocum/L0-timeseries/') + '/'
profiledir = str(example_dir / 'example-slocum/L0-profiles/')

outname = slocum.binary_to_timeseries(
    binarydir,
    cacdir,
    tsdir,
    deploymentyaml,
    search='*.[s|t]bd',
    profile_filt_time=20,
    profile_min_time=20,
)
ncprocess.extract_timeseries_profiles(outname, profiledir, deploymentyaml, force=True)

output = xr.open_dataset(outname)

# ---------------------------------------------------------------------------
# Load expected YAML summary
# ---------------------------------------------------------------------------
expected_yaml_path = (
    library_dir
    / 'tests/expected/example-slocum/L0-timeseries/dfo-rosie713-20190615.yml'
)
with open(expected_yaml_path) as f:
    expected = yaml.safe_load(f)




def _run_compliance(path, checker_names):
    """Return the compliance_checker JSON result dict for *path*."""
    check_suite = CheckSuite()
    check_suite.load_all_available_checkers()
    report = example_dir / 'report.json'
    ComplianceChecker.run_checker(
        str(path),
        checker_names,
        verbose=0,
        criteria='normal',
        output_filename=str(report),
        output_format='json',
    )
    with open(report) as fp:
        return json.load(fp)


# ---------------------------------------------------------------------------
# Timeseries tests
# ---------------------------------------------------------------------------

def test_variable_list():
    assert sorted(output.variables) == sorted(expected['variables'])


def test_global_attrs():
    out_attrs = {k: str(v) for k, v in output.attrs.items() if k not in SKIP_ATTRS}
    exp_attrs = {k: v for k, v in expected['attrs'].items() if k not in SKIP_ATTRS}
    assert out_attrs == exp_attrs


@pytest.mark.parametrize('var', sorted(expected['variables']))
def test_variable_attrs(var):
    exp_attrs = expected['variables'][var]['attrs']
    out_attrs = {k: str(v) for k, v in output[var].attrs.items()}
    assert out_attrs == exp_attrs


@pytest.mark.parametrize('var', sorted(expected['variables']))
def test_variablestats(var):
    vals = to_float(output[var])
    actual = stats(vals)
    exp = expected['variables'][var]

    assert actual['n_valid'] == exp['n_valid'], f'{var}: n_valid mismatch'
    assert actual['n_nan'] == exp['n_nan'], f'{var}: n_nan mismatch'

    for stat in ('min', 'max', 'mean', 'std'):
        if stat not in exp:
            continue
        np.testing.assert_allclose(
            actual[stat], exp[stat], rtol=RTOL, atol=ATOL,
            err_msg=f'{var}: {stat} mismatch',
        )


@pytest.mark.parametrize('var', sorted(expected['variables']))
def test_variable_diffstats(var):
    """Time-derivative statistics catch temporal scrambling."""
    vals = to_float(output[var])
    actual = stats(vals)
    exp = expected['variables'][var]

    for stat in ('diff_mean', 'diff_std'):
        if stat not in exp:
            continue
        np.testing.assert_allclose(
            actual[stat], exp[stat], rtol=RTOL, atol=ATOL,
            err_msg=f'{var}: {stat} mismatch',
        )


def test_time_monotonic():
    """Time must be strictly increasing."""
    t = output['time'].values.astype('datetime64[ns]').astype('float64')
    assert np.all(np.diff(t) > 0), 'time is not monotonically increasing'


# ---------------------------------------------------------------------------
# Compliance tests (ported from test_slocum.py)
# ---------------------------------------------------------------------------

# @pytest.mark.xfail(reason='profiles not fully compliant')
def _compliance_msgs(cc_result, priority):
    """Extract messages from all checks at a given priority level."""
    key = {3: 'high_priorities', 2: 'medium_priorities', 1: 'low_priorities'}[priority]
    return [
        f"{check['name']}: {msg}"
        for check in cc_result.get(key, [])
        for msg in check.get('msgs', [])
    ]


def test_profiles_compliant():
    path = Path(profiledir) / 'dfo-rosie713-20190620T1313.nc'
    cc_data = _run_compliance(path, ['gliderdac', 'cf:1.8'])
    for checker in ['gliderdac', 'cf:1.8']:
        result = cc_data[checker]
        assert result['high_count'] == 0, \
            f"{checker} high priority errors:\n" + "\n".join(_compliance_msgs(result, 3))
        assert result['medium_count'] == 0, \
            f"{checker} medium priority errors:\n" + "\n".join(_compliance_msgs(result, 2))
        assert result['low_count'] == 0, \
            f"{checker} low priority errors:\n" + "\n".join(_compliance_msgs(result, 1))


#@pytest.mark.xfail(strict=False,
#                   reason='compliance_checker result varies across versions')
def test_timeseries_compliant():
    cc_data = _run_compliance(outname, ['cf:1.8'])
    result = cc_data['cf:1.8']
    assert result['high_count'] == 0, \
        "cf:1.8 high priority errors:\n" + "\n".join(_compliance_msgs(result, 3))
    assert result['medium_count'] == 0, \
        "cf:1.8 medium priority errors:\n" + "\n".join(_compliance_msgs(result, 2))
    assert result['low_count'] == 0, \
        "cf:1.8 low priority errors:\n" + "\n".join(_compliance_msgs(result, 1))
