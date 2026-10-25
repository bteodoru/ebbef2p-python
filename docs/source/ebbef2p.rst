ebbef2p package
===============

The names most users need are available directly from the package:
:class:`~ebbef2p.model.Model`, :func:`~ebbef2p.solver.solve`,
:class:`~ebbef2p.vlasov.VlasovSoil` and :func:`~ebbef2p.vlasov.solve_vlasov`.

Model
-----

.. automodule:: ebbef2p.model
   :members:

Solver
------

.. automodule:: ebbef2p.solver
   :members:

Results
-------

.. automodule:: ebbef2p.results
   :members:

Modified Vlasov model
---------------------

.. automodule:: ebbef2p.vlasov
   :members:

Mesh
----

.. automodule:: ebbef2p.mesh
   :members:

Element matrices
----------------

.. automodule:: ebbef2p.elements
   :members:

Exceptions
----------

All errors raised by the package derive from
:class:`~ebbef2p.exceptions.EBBEF2pError` and also from the matching
built-in exception, so they can be caught either way.

.. automodule:: ebbef2p.exceptions
   :members:
   :show-inheritance:

Plotting
--------

.. automodule:: ebbef2p.plotting
   :members:
