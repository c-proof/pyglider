from setuptools import setup, find_packages

# Get the version from versioneer
__version__ = "0.0.4"

setup(name="pyglider",
      version=__version__,
      description="Glider data to netCDF translation in python",
      author="Jody Klymak",
      author_email="jklymak@gmail.com",
      url="https://pyglider.readthedocs.io",
      packages=find_packages(exclude=['tests']),
      python_requires='>=3.6',
      install_requires=[
        "xarray",
        "dask",
        "netcdf4",
        "gsw",
        "scipy",
        "bitstring",
        "pooch"
      ],
      license='Apache',
      extras_require={
        "code_style": ["flake8<3.8.0,>=3.7.0", "black", "pre-commit==1.17.0"],
        "testing": ["pytest"],
        "docs": ["pydata-sphinx-theme", "numpydoc", "autoapi", "myst-parser",
                 "sphinx"]
      },
      zip_safe=True
      )
