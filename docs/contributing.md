# Contributing

## Create a fork

```{note}
TODO: Create contributing page noting:
- the forking process (trimmed down)
```

## Development environment

Assuming that you have Anaconda/Miniconda installed, you can create a new environment for development:

```bash
conda create -n pyglider-dev python=3.9
conda activate pyglider-dev
```

Then install the dependencies:

```bash
conda env update -f environment.yml -n pyglider-dev
```

And install the package in editable mode:

```bash
pip install -e . --no-build-isolation --no-deps
```

Done!

---

From here, you can make the changes you want, and add tests. When you are ready, you can create a pull request into the codebase.

## Running tests

Once the development environment is set up, you can run the tests with:

```bash
pytest
```

## Building documentation

```{note}
TODO
```

## [Optional] Running pre-commit

We use pre-commit [TODO link] to run tooling on the code to make sure that it adheres to standards that we have adopted in the codebase. This is done automatically in the cloud on all pull requests, however this can also be done locally.

Pre-commit is already a dependency in the `environment.yml` file, so you only have to do `pre-commit install` in the repository (this will install the hooks in the repository so that they run when you commit changes).

That's it! If you want to, you can manually run the hooks by doing `pre-commit run --all-files` (or `pre-commit run` if you only want to run on staged files).
