# Development

## Installation for development

PyGlider uses [pixi](https://pixi.sh/) to manage environments. To get started, clone the repo and install in editable mode:

```bash
git clone https://github.com/c-proof/pyglider.git
cd pyglider
pixi install
```

This installs pyglider itself as an editable package (changes to source files take effect immediately without reinstalling). The default environment includes the library and its dependencies. To also get the test dependencies:

```bash
pixi install -e test
```

If you prefer conda/pip without pixi, install the dependencies with conda and then do an editable pip install:

```bash
conda create -n pyglider-dev
conda activate pyglider-dev
conda install -c conda-forge dask netcdf4 xarray numpy scipy gsw
pip install -e .
```

## Running the tests

PyGlider uses [pytest](https://docs.pytest.org/). To run the tests with pixi:

```bash
pixi run -e test pytest tests/
```

Or, if your environment already has the dependencies installed:

```bash
pytest tests/
```

Tests run the full processing pipeline on example data and compare the output against stored NetCDF golden files in `tests/expected/`. A test failure means the output has changed — either a regression or an intentional change that requires updating the golden files (see below).

## How the regression tests work

Each golden file pair in `tests/expected/` consists of:

- A cleaned `.nc` file — the reference dataset used for comparison.  Dynamic
  global attributes (`date_created`, `date_issued`, `history`) are stripped so
  the file is stable across runs.
- A `.cdl` header file (`ncdump -h` output) alongside it for human-readable
  diffs in git.

Tests compare fresh pipeline output against the golden `.nc` files using
`xarray.testing.assert_identical` (via `nc_test_helpers.assert_datasets_equal`),
stripping dynamic attributes from the actual output before comparison.

The golden files live under `tests/expected/`, organized by glider type and pipeline stage:

```
tests/expected/
    example-seaexplorer/
        L0-timeseries/
            dfo-eva035-20190718.nc
            dfo-eva035-20190718.cdl
            dfo-eva035-20190718_adjusted.nc
            dfo-eva035-20190718_adjusted.cdl
        L0-gridfiles/
            dfo-eva035-20190718_grid_adjusted.nc
            dfo-eva035-20190718_grid_adjusted.cdl
    example-seaexplorer-raw/
        L0-timeseries/
            dfo-bb046-20200908.nc
            dfo-bb046-20200908.cdl
    example-slocum/
        L0-timeseries/
            dfo-rosie713-20190615.nc
            dfo-rosie713-20190615.cdl
            dfo-rosie713-20190615_adjusted.nc
            dfo-rosie713-20190615_adjusted.cdl
        L0-gridfiles/
            dfo-rosie713-20190615_grid_adjusted.nc
            dfo-rosie713-20190615_grid_adjusted.cdl
```

The test files (`tests/test_*_nc.py`) each run the pipeline once at module
load time, writing fresh output to `tests/example-data/example-*/L0-timeseries-test/`
(and equivalent `-test` directories for other pipeline stages), then compare
that output against the golden NC files in `tests/expected/`.  The
`process_adjusted` tests are an exception — they write directly into
`tests/expected/` and clean the files in place.

## Updating golden files after an intentional change

If you make a change to pyglider that intentionally alters the output — new variable, changed attribute, fixed a calculation — the golden files need to be updated.

**Step 1: Regenerate the golden files:**

```bash
python tests/_generate_expected_cdl.py
```

This re-runs the processing pipelines, writes cleaned NC golden files to
`tests/expected/`, and generates `.cdl` headers alongside each one.  For the
slocum L0 timeseries pipeline (which requires `dbdreader`), the golden file is
read from a pre-existing NC in `tests/example-data/example-slocum/L0-timeseries/`;
run `tests/example-data/example-slocum/process_deploymentRealTime.py` to
refresh that file if needed.

**Step 2: Re-run the tests** to confirm they now pass:

```bash
pixi run -e test pytest tests/
```

**Step 3: Inspect the diff** before committing:

```bash
git diff tests/expected/
```

The diff should show only the changes you intended. If unrelated variables or datasets changed, investigate before committing.

## Test file inventory

| Test file | What it covers |
|-----------|----------------|
| `test_seaexplorer_nc.py` | SeaExplorer NRT sub and raw delayed L0 timeseries; interpolation behaviour |
| `test_slocum_nc.py` | Slocum L0 timeseries and profiles; CF/GliderDAC compliance |
| `test_seaexplorer_og10_nc.py` | SeaExplorer OG 1.0 timeseries; OG/CF compliance |
| `test_slocum_og10_nc.py` | Slocum OG 1.0 timeseries; OG/CF compliance |
| `test_process_adjusted_seaexplorer_nc.py` | SeaExplorer `process_adjusted` timeseries and gridfiles |
| `test_process_adjusted_slocum_nc.py` | Slocum `process_adjusted` timeseries and gridfiles |
| `test_seaexplorer.py` | SeaExplorer-specific unit tests |
| `test_utils.py` | Utility function unit tests |

## Making a release

Releases are published to PyPI automatically when a version tag is pushed. The
GitHub Actions workflow builds the package, publishes to TestPyPI, then PyPI,
and creates a GitHub Release with auto-generated notes.

The version number is derived from git tags via `setuptools-scm` — there is no
version to manually bump.

**Step 1: Tag the commit you want to release:**

```bash
git tag v1.2.3
git push --tags
```

That's it. The workflow triggers on any `v*` tag.

**TestPyPI** runs first as a smoke test. If it fails (e.g. the package already
exists at that version on TestPyPI), you can re-tag after fixing, but note that
TestPyPI does not allow re-uploading the same version — use a new tag.

**Trusted Publishing** must be configured on both PyPI and TestPyPI for the
`pypi` and `testpypi` GitHub Environments respectively. If the workflow's
publish step fails with a permissions error, check that the environment is
configured in the repo settings and that the PyPI trusted publisher entry
matches the repo and workflow file name.

## Numerical tolerances

Comparisons use `xarray.testing.assert_identical`, which requires exact bit-for-bit equality of all variable values and attributes. Time monotonicity is checked as a separate test. The `date_created`, `date_issued`, and `history` global attributes are excluded from comparison because they change on every run.
