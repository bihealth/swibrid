# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "swibrid"
copyright = "2024, Benedikt Obermayer"
author = "Benedikt Obermayer"
release = "2024"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinxcontrib.autoprogram", "sphinx.ext.autodoc"]

templates_path = ["_templates"]
exclude_patterns = []

import os
import sys

# sys.path.insert(0, os.path.abspath(os.path.join('..', '..','..', 'swibrid')))
sys.path.insert(0, os.path.abspath(os.path.join("..", "..", "swibrid")))
sys.path.insert(
    0,
    os.path.abspath(
        os.path.join(
            "..",
            "..",
        )
    ),
)
sys.path.insert(0, os.path.abspath(os.path.join("..")))
# sys.path.insert(0, os.path.abspath(os.path.join('..', 'swibrid')))

autodoc_member_order = "bysource"


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
