Getting started
===============

A problem is described by a :class:`~ebbef2p.model.Model`, solved with
:func:`~ebbef2p.solver.solve`, and inspected through the returned
:class:`~ebbef2p.results.Solution`.

Defining a model
----------------

.. code-block:: python

    from ebbef2p import Model, solve

    model = Model()
    model.add_beam(0, 10, E=30e6, I=0.01)          # one or more segments
    model.add_foundation(0, 10, k=5000, g=2000)     # Winkler k, shear g = 2t
    model.add_point_load(5, force=100)              # force and/or moment
    model.add_distributed_load(0, 10, 20)           # uniform; (x0, x1, q0, q1) for linear

All ``add_*`` methods return the model, so calls can be chained.

Beam segments
    :meth:`~ebbef2p.model.Model.add_beam` adds a segment with constant ``E``
    and ``I``. Segments must cover the beam without gaps or overlaps; this
    is how a beam of stepped section is modelled.

Foundation
    :meth:`~ebbef2p.model.Model.add_foundation` adds a segment on which the
    first parameter ``k`` and the second parameter ``g`` vary linearly. Each
    may be a scalar or a ``(start, end)`` pair. Overlapping segments add up,
    so ``k`` and ``g`` can be described independently. ``g`` is 2t in
    Vlasov's notation and G\ :sub:`p` in Pasternak's; with ``g = 0`` the
    foundation is of Winkler type.

Soil beyond the ends
    ``Model(soil_beyond_ends=True)`` lets the foundation continue
    indefinitely beyond both beam ends. Its effect is an equivalent spring
    :math:`\sqrt{k\,g}` at each end (see :doc:`../theory/fe_modeling`).

Supports
    :meth:`~ebbef2p.model.Model.add_support` prescribes the displacement
    and/or the rotation at a point. The default is a simple support
    (``w = 0``, rotation free); ``theta=0`` clamps the beam; a non-zero
    ``w`` prescribes a settlement.

Solving and results
-------------------

.. code-block:: python

    solution = solve(model, n_elements=100)

    solution.x, solution.w, solution.theta   # nodal coordinates and values
    solution.deflection([2.5, 7.5])           # anywhere along the beam

    M = solution.bending_moment()             # values at element ends
    V = solution.shear_force()
    M.max(), M.min(), M.at(5, side='left')

    p = solution.reactive_pressure()
    p.x, p.values          # distributed soil pressure
    p.concentrated         # [(x, force), ...] at the ends and at jumps of g
    p.resultant            # total soil reaction

    solution.support_reactions()              # [(x, force, moment), ...]
    solution.equilibrium()                    # applied vs reactions

Internal forces and pressures are returned at both ends of every element,
so jumps at nodes (under a point load, where the section or the foundation
changes) are preserved. :meth:`~ebbef2p.results.ElementEndValues.at` picks
the value on either side of a node.

Modified Vlasov foundation
--------------------------

When the soil is described as an elastic layer rather than by ``k`` and
``g``, use :func:`~ebbef2p.vlasov.solve_vlasov`:

.. code-block:: python

    from ebbef2p import VlasovSoil, solve_vlasov

    soil = VlasovSoil(E_top=144, nu=0.2, depth=25)    # optional E_bottom
    result = solve_vlasov(model_without_foundation, soil)
    result.gamma, result.k, result.g, result.solution

Quick plots
-----------

.. code-block:: python

    from ebbef2p.plotting import plot_solution
    plot_solution(solution)
