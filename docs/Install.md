# Installation

## conda/pip

PyGlider depends on `dask` and `netcdf4`, both of which can be tricky to install using `pip`,
hence we recommend these be installed with [`conda`](https://www.anaconda.com/). To install
PyGlider, create an environment, and do

```
conda create -n gliderwork
conda activate gliderwork
conda install -c conda-forge pyglider
```

## Editable installation

Follow the instructions in our [contributing page](./contributing.md).
