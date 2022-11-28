import os
import sys
sys.path.insert(0, os.path.abspath('../rbatools'))
# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'rbatools'
copyright = '2022 INRAE - MaIAGE - France'
author = 'Bodeit, O., Liebermeister, W. and Goelzer A.'

release = '1.0'
version = '1.0.0'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.napoleon'
]

napoleon_google_docstring = False
napoleon_numpy_docstring = True


intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'
