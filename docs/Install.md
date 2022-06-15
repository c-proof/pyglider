# Installation

## Dependencies

- install conda/miniconda
- clone the repository onto your computer and cd into the directory
- from inside the cloned repository directory run
`conda env create -f environment.yml`, it will read environment.yml and
automatically install the dependencies, and create environment `pyglider`,
(though you can change the name by editing `environment.yml`)

## Editable installation

If you want to be able to edit the files in `pyglider/pyglider` then install
as above, and then do `pip install -e .`.  That will re-install pyglider
with links to the local directory, so you can edit the library files.
If you do so consider making a pull-request with your changes!