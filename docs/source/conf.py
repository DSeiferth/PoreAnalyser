# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'PoreFinding'
copyright = '2023, Seiferth'
author = 'David Seiferth'

release = '0.1'
version = '0.1.0'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
]

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

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

#import PoreFinding
#Traceback (most recent call last):
#  File "/home/docs/checkouts/readthedocs.org/user_builds/porefinding/envs/latest/lib/python3.10/site-packages/sphinx/config.py", line 356, in eval_config_file
#    exec(code, namespace)  # NoQA: S102
#  File "/home/docs/checkouts/readthedocs.org/user_builds/porefinding/checkouts/latest/docs/source/conf.py", line 41, in <module>
#    import PoreFinding
#  File "/home/docs/checkouts/readthedocs.org/user_builds/porefinding/checkouts/latest/PoreFinding/__init__.py", line 8, in <module>
#    from .porefinding import PoreAnalysis
#  File "/home/docs/checkouts/readthedocs.org/user_builds/porefinding/checkouts/latest/PoreFinding/porefinding.py", line 6, in <module>
#    import hole_analysis as hole_analysis
#ModuleNotFoundError: No module named 'hole_analysis'

