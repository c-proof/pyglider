# Testing

## Running the tests

PyGlider uses [pytest](https://docs.pytest.org/). The test environment is managed with [pixi](https://pixi.sh/):

```bash
pixi run -e test pytest tests/
```

Or, if your environment already has the dependencies installed:

```bash
pytest tests/
```

Tests run the full processing pipeline on example data and compare the output against stored YAML "golden files" in `tests/expected/`. A test failure means the output has changed — either a regression or an intentional change that requires updating the golden files (see below).

## How the YAML regression tests work

Rather than committing binary NetCDF golden files to the repo, the tests summarize each output dataset into a compact YAML file and compare that instead. Each YAML file captures:

- the list of variables
- global attributes (excluding dynamic fields like `date_created` and `history`)
- per-variable CF attributes
- per-variable statistics: `min`, `max`, `mean`, `std`, `n_valid`, `n_nan`
- per-variable time-derivative statistics: `diff_mean`, `diff_std`

This makes diffs human-readable: if a processing change shifts the mean salinity by 0.01 PSU, you see exactly that in the YAML diff rather than a binary blob change.

The golden files live under `tests/expected/`, organized by glider type and pipeline stage:

```
tests/expected/
    example-seaexplorer/
        L0-timeseries/
            dfo-eva035-20190718.yml
            dfo-eva035-20190718_adjusted.yml
        L0-gridfiles/
            dfo-eva035-20190718_grid_adjusted.yml
    example-seaexplorer-raw/
        L0-timeseries/
            dfo-bb046-20200908.yml
    example-slocum/
        L0-timeseries/
            dfo-rosie713-20190615.yml
            dfo-rosie713-20190615_adjusted.yml
        L0-gridfiles/
            dfo-rosie713-20190615_grid_adjusted.yml
```

The test files (`tests/test_*_yaml.py`) each run the pipeline once at module load time, then use parametrized pytest tests to check each variable and attribute.

## Updating golden files after an intentional change

If you make a change to pyglider that intentionally alters the output — new variable, changed attribute, fixed a calculation — the YAML golden files need to be updated to reflect the new expected values.

**Step 1: Run the tests** to generate up-to-date output NC files:

```bash
pixi run -e test pytest tests/
```

The tests will fail (that's expected), but the pipeline will have written fresh NC files to `tests/example-data/`.

**Step 2: Regenerate the YAML golden files** from those NC files:

```bash
python tests/_generate_expected_yaml.py
```

This reads each NC file, computes the same summary statistics used by the tests, and overwrites the YAML files in `tests/expected/`.

**Step 3: Re-run the tests** to confirm they now pass:

```bash
pixi run -e test pytest tests/
```

**Step 4: Inspect the diff** before committing:

```bash
git diff tests/expected/
```

The diff should show only the changes you intended. If unrelated variables or datasets changed, investigate before committing.

## Test file inventory

| Test file | What it covers |
|-----------|----------------|
| `test_seaexplorer_yaml.py` | SeaExplorer NRT sub and raw delayed L0 timeseries; interpolation behaviour |
| `test_slocum_yaml.py` | Slocum L0 timeseries and profiles; CF/GliderDAC compliance |
| `test_process_adjusted_seaexplorer_yaml.py` | SeaExplorer `process_adjusted` timeseries and gridfiles |
| `test_process_adjusted_slocum_yaml.py` | Slocum `process_adjusted` timeseries and gridfiles |
| `test_seaexplorer.py` | SeaExplorer-specific unit tests |
| `test_utils.py` | Utility function unit tests |

## Numerical tolerances

Statistics comparisons use relative tolerance `rtol=1e-5` and `atol=0.0`. Integer counts (`n_valid`, `n_nan`) are compared exactly. Time monotonicity is checked as a separate test. The `date_created`, `date_issued`, and `history` global attributes are excluded from comparison because they change on every run.
