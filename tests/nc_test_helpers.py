"""
Shared helpers and constants for CDL/NC-based regression tests.

Imported by test_slocum_yaml.py, test_seaexplorer_yaml.py,
test_process_adjusted_slocum_yaml.py,
test_process_adjusted_seaexplorer_yaml.py, and the OG 1.0 variants.
"""

import json
from pathlib import Path

import xarray as xr
from compliance_checker.runner import CheckSuite, ComplianceChecker

# Global attributes excluded from comparison (dynamic / version-stamped fields)
SKIP_ATTRS = {'date_created', 'date_issued', 'history'}


def assert_datasets_equal(actual: xr.Dataset, expected: xr.Dataset) -> None:
    """Compare two datasets, ignoring dynamic global attributes.

    Strips SKIP_ATTRS from *actual* before calling xarray.testing.assert_identical.
    Golden NC files are already stored without those attrs (see _generate_expected_cdl.py).
    """
    a = actual.copy()
    for attr in SKIP_ATTRS:
        a.attrs.pop(attr, None)
    xr.testing.assert_identical(a, expected)


def run_compliance(path, checker_names):
    """Return the compliance_checker JSON result dict for *path*."""
    check_suite = CheckSuite()
    check_suite.load_all_available_checkers()
    report = Path(path).parent / 'report.json'
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


def compliance_msgs(cc_result, priority):
    """Extract messages from all checks at a given priority level."""
    key = {3: 'high_priorities', 2: 'medium_priorities', 1: 'low_priorities'}[priority]
    return [
        f"{check['name']}: {msg}"
        for check in cc_result.get(key, [])
        for msg in check.get('msgs', [])
    ]
