Installation
============

|name| is available on the |pypi_long| (PyPI) and can be installed with the
Python package installer pip_:

.. code-block:: console

    pip install ebbef2p-python

To install the development version from GitHub:

.. code-block:: console

    pip install git+https://github.com/bteodoru/ebbef2p-python.git

|name| requires Python 3.12 or newer and depends on NumPy and SciPy.
Matplotlib is needed only for :func:`ebbef2p.plotting.plot_solution`; it is
installed with the ``plot`` extra:

.. code-block:: console

    pip install "ebbef2p-python[plot]"

Development setup
-----------------

The project is managed with uv_. After cloning the repository,

.. code-block:: console

    uv sync                 # creates .venv with the package and dev tools
    uv run pytest           # unit and verification tests
    uv run ruff check .
    uv run sphinx-build -W -b html docs/source docs/build

The tests are split into ``tests/unit`` (API behaviour, validation,
internal consistency) and ``tests/verification`` (comparisons with
analytical solutions, published results, numerical integration and global
equilibrium).

.. _uv: https://docs.astral.sh/uv/

.. seealso::

    How to get started with pip_: https://pip.pypa.io/en/stable/quickstart/
