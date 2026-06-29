"""
Shared helpers and constants for CDL/NC-based regression tests.

Imported by test_seaexplorer_nc.py, test_slocum_nc.py,
test_seaexplorer_og10_nc.py, test_slocum_og10_nc.py,
test_process_adjusted_seaexplorer_nc.py, and test_process_adjusted_slocum_nc.py.
"""

import json
from pathlib import Path

import xarray as xr
from compliance_checker.runner import CheckSuite, ComplianceChecker


def assert_datasets_equal(actual: xr.Dataset, expected: xr.Dataset) -> None:
    """Compare two datasets on data values only, ignoring attributes.

    Attributes are tested via the .cdl golden files (git-diffable).
    Uses assert_allclose to tolerate minor cross-platform float differences
    (e.g. from different GSW/scipy versions between local and CI).
    """
    xr.testing.assert_allclose(actual, expected)


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
