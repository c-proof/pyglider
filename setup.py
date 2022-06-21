from setuptools import setup, find_packages

# Get the version from versioneer
__version__ = "0.0.2"

with open("./docs/requirements.txt", "r") as handle:
    requirements_doc = [ld.strip() for ld in handle.read().splitlines() if ld.strip()]

with open("./tests/requirements.txt", "r") as handle:
    requirements_dev = [ld.strip() for ld in handle.read().splitlines() if ld.strip()]

with open("./requirements.txt", "r") as handle:
    requirements = [ld.strip() for ld in handle.read().splitlines() if ld.strip()]

setup(name="pyglider",
      version=__version__,
      description="Glider data to netCDF translation in python",
      author="Jody Klymak",
      author_email="jklymak@gmail.com",
      url="https://pyglider.readthedocs.io",
      packages=find_packages(exclude=['tests']),
      python_requires='>=3.6',
      install_requires=requirements,
      extras_require={
        "code_style": ["flake8<3.8.0,>=3.7.0", "black", "pre-commit==1.17.0"],
        "testing": requirements_dev,
        "docs": requirements_doc,
      },
      zip_safe=True
      )
