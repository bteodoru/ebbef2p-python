"""Sphinx configuration for EBBEF2p."""

from importlib.metadata import version as _version

NAME = 'EBBEF2p'
NAME_LONG = 'Euler-Bernoulli Beams on Two-Parameter Elastic Foundations'
NAME_LC = NAME.lower()
NAME_CLI = NAME_LC
AUTHOR = 'Bogdan Teodoru'
LICENCE_NAME = "MIT"

project = 'ebbef2p-python'
copyright = '2020-2026, Bogdan Teodoru'
author = AUTHOR
release = _version('ebbef2p-python')
version = '.'.join(release.split('.')[:2])

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.coverage',
    'sphinx.ext.napoleon',
    'sphinx.ext.doctest',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
]
autodoc_member_order = 'bysource'
master_doc = 'index'
language = 'en'
exclude_patterns = []

rst_prolog = f"""
.. |name| replace:: {NAME}
.. |name_long| replace:: {NAME_LONG}
.. |name_cli| replace:: ``{NAME_CLI}``
.. |eb| replace:: Euler-Bernoulli
.. |name_bold| replace:: **{NAME}**
.. |dof| replace:: d.o.f.
.. |dof_long| replace:: degrees of freedom
.. _FEM: https://en.wikipedia.org/wiki/Finite_element_method
.. _FE: https://en.wikipedia.org/wiki/Finite_element_method
.. _PyPI: https://pypi.org/project/ebbef2p-python/
.. |pypi_long| replace:: Python Package Index
.. _pip: https://pip.pypa.io/en/stable/
.. |author1| replace:: {AUTHOR}
.. |license| replace:: {LICENCE_NAME}
"""

html_theme = 'alabaster'
html_theme_options = {
    'show_powered_by': False,
    'github_user': 'bteodoru',
    'github_repo': 'ebbef2p-python',
    'github_banner': False,
    'show_related': True,
    'note_bg': '#FFF59C',
}
