"""Utilities specific to the test suite."""

from pathlib import Path

LIBRARY_DIR = Path(__file__).parent.parent.absolute()
EXAMPLE_DIR = LIBRARY_DIR / 'tests/example-data/'
