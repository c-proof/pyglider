"""
Shared helpers and constants for YAML-based regression tests.

Imported by test_slocum_yaml.py, test_seaexplorer_yaml.py,
test_process_adjusted_slocum_yaml.py,
test_process_adjusted_seaexplorer_yaml.py, and the OG 1.0 variants.
"""

import json
from pathlib import Path

import numpy as np
from compliance_checker.runner import CheckSuite, ComplianceChecker

# Tolerances for per-variable statistic comparisons
RTOL = 1e-8
ATOL = 0.0

# Global attributes excluded from comparison (dynamic / version-stamped fields)
SKIP_ATTRS = {'date_created', 'date_issued', 'history'}


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


def to_float(da):
    """Return a flat float64 array, converting datetime64 to ns-epoch float.
    Returns an empty array for non-numeric, non-datetime variables (e.g. strings)."""
    vals = da.values
    if not np.issubdtype(vals.dtype, np.number) and not np.issubdtype(vals.dtype, np.datetime64):
        return np.array([], dtype='float64')
    if np.issubdtype(vals.dtype, np.datetime64):
        vals = vals.astype('datetime64[ns]').astype('float64')
    return vals.flatten().astype('float64')


def stats(vals):
    """Compute summary statistics dict for a flat float64 array."""
    valid = vals[~np.isnan(vals)]
    diffs = np.diff(valid)
    s = {
        'n_valid': int(len(valid)),
        'n_nan': int(np.sum(np.isnan(vals))),
    }
    if len(valid):
        s.update(
            min=float(np.min(valid)),
            max=float(np.max(valid)),
            mean=float(np.mean(valid)),
            std=float(np.std(valid)),
        )
    if len(diffs):
        s.update(
            diff_mean=float(np.mean(diffs)),
            diff_std=float(np.std(diffs)),
        )
    return s
