# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

from pyglider._version import __version__

# -- Project information -----------------------------------------------------

project = 'PyGlider'
copyright = '2022-, PyGlider team'
author = 'PyGlider team'

# The full version, including alpha/beta/rc tags
release = __version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "numpydoc",
    "myst_parser",
    'sphinx.ext.autodoc',
    'sphinx.ext.inheritance_diagram',
    'autoapi.sphinx',
]

extensions.append('sphinx.ext.intersphinx')

intersphinx_mapping = {
  'xarray': ('http://xarray.pydata.org/en/stable/', None),
  'python': ('https://docs.python.org/3/', None),
  }

autoapi_modules = {'pyglider': None}

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "pydata_sphinx_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_logo = "_static/PyGliderHorizontal.svg"

html_context = {
    # "github_url": "https://github.com", # or your GitHub Enterprise interprise
    "github_user": "c-proof",
    "github_repo": "pyglider",
    "doc_path": "docs/",
}
html_theme_options = {
    "icon_links": [
        {
            # Label for this link
            "name": "GitHub",
            # URL where the link will redirect
            "url": "https://github.com/c-proof/pyglider",  # required
            # Icon class (if "type": "fontawesome"), or path to local image (if "type": "local")
            "icon": "fab fa-github-square",
            # The type of image to be used (see below for details)
            "type": "fontawesome",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/pyglider",
            "icon": "fas fa-box",
        },
        # {   "name": "conda-forge",
        #     "url": "https://anaconda.org/conda-forge/pyglider",
        #     "icon": "fas fa-box"
        # }
   ]}
