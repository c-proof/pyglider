import re
from setuptools import setup, find_packages

# Get the version from versioneer
__version__ = "0.0"

setup(name="pyglider",
      version=__version__,
      description="Glider processing toolbox in python",
      author="jklymak",
      author_email="jklymak@gmail.com",
      url="https://github.com/jklymak/pyglider.git",
      packages=find_packages(exclude=['tests']),
      python_requires='>=3.6',
      install_requires=[
          "matplotlib",
          "numpy",
          "xarray",
          "dask",
          "bitstring",
          "netcdf4",
          "geojson",
          "seawater",
          "scipy",
      ],
      zip_safe=True
      )
