"""
Shared helpers and constants for YAML-based regression tests.

Imported by test_slocum_yaml.py, test_seaexplorer_yaml.py,
test_process_adjusted_slocum_yaml.py, and
test_process_adjusted_seaexplorer_yaml.py.
"""

import numpy as np

# Tolerances for per-variable statistic comparisons
RTOL = 1e-8
ATOL = 0.0

# Global attributes excluded from comparison (dynamic / version-stamped fields)
SKIP_ATTRS = {'date_created', 'date_issued', 'history'}


def to_float(da):
    """Return a flat float64 array, converting datetime64 to ns-epoch float."""
    vals = da.values
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
