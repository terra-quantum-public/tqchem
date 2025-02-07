# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sys
from pathlib import Path

import tomllib

tqchem_root = Path(__file__).parents[2].resolve()
sys.path.insert(0, tqchem_root)

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

with open(tqchem_root / "pyproject.toml", "rb") as f:
    pyproject_data = tomllib.load(f)

project = pyproject_data["project"]["name"]
release = pyproject_data["project"]["version"]
author = ", ".join([author["name"] for author in pyproject_data["project"]["authors"]])
copyright = f"2024, {author}"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "autoapi.extension",
    "nbsphinx",
    "nbsphinx_link",
    "sphinx_rtd_theme",
]
templates_path = ["_templates"]
exclude_patterns = ["_build", "**.ipynb_checkpoints"]

autoapi_type = "python"
autoapi_dirs = ["../../tqchem/"]
autoapi_options = [  # default without private functions
    "members",
    "undoc-members",
    "show-inheritance",
    "show-module-summary",
    "special-members",
    "imported-members",
]
autoapi_ignore = ["*/__main__.py"]
autoapi_file_patterns = ["*.pyi", "*.py"]

html_theme = "sphinx_rtd_theme"
