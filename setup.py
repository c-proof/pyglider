from setuptools import setup, find_packages

# Get the version from versioneer
__version__ = "0.0.2"

setup(name="pyglider",
      version=__version__,
      description="Glider data to netCDF translation in python",
      author="Jody Klymak",
      author_email="jklymak@gmail.com",
      url="https://github.com/c-proof/pyglider.git",
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
          "gsw",
          "scipy",
          "pooch",
          "pyyaml"
      ],
      zip_safe=True
      )
