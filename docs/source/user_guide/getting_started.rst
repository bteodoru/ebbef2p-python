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

    For a rectangular section, :func:`~ebbef2p.sections.rectangular_section`
    gives ``I = b*h**3/12``. Take ``b`` as the width in contact with the
    soil; for any other section enter the equivalent height
    ``h_eq = (12*I/b)**(1/3)`` computed with that width. The same ``b``
    turns soil parameters given per unit area into the per-unit-length
    values of the model, ``k = k_s*b``, so a width chosen only to reproduce
    ``I`` would give a wrong ``k``.

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
the value on either side of a node, and
:meth:`~ebbef2p.results.ElementEndValues.interpolate` gives values anywhere
on the beam, linear between the ends of each element (``M.interpolate(2.3)``,
``M.interpolate(5, side='right')``). Between nodes this is an interpolation,
so its accuracy depends on the mesh; for the reactive pressure it covers the
distributed part only.

Modified Vlasov foundation
--------------------------

When the soil is described as an elastic layer rather than by ``k`` and
``g``, use :func:`~ebbef2p.vlasov.solve_vlasov`:

.. code-block:: python

    from ebbef2p import VlasovSoil, solve_vlasov

    soil = VlasovSoil(E_top=144, nu=0.2, depth=25)    # optional E_bottom
    result = solve_vlasov(model_without_foundation, soil)
    result.gamma, result.k, result.g, result.solution

The layer continues beyond the beam ends only if the model was created with
``Model(soil_beyond_ends=True)``.

Saving and loading models
-------------------------

:mod:`ebbef2p.io` turns a model into a JSON-ready dictionary and back:

.. code-block:: python

    import json
    from ebbef2p.io import model_from_dict, model_to_dict

    with open("beam.json", "w") as f:
        json.dump(model_to_dict(model), f, indent=2)

    with open("beam.json") as f:
        model = model_from_dict(json.load(f))

The model is rebuilt through the ``add_*`` methods and checked as a whole,
so errors are those of building it in code, with the path of the offending
field (``beams[1].E: expected a number, got '2000'``); unknown keys are
errors. ``model_from_dict(data, validate=False)`` reloads an incomplete
model. :func:`~ebbef2p.io.vlasov_soil_to_dict` and
:func:`~ebbef2p.io.vlasov_soil_from_dict` do the same for a
:class:`~ebbef2p.vlasov.VlasovSoil`, and
:func:`~ebbef2p.io.model_json_schema` gives a JSON Schema of the format in
which every field is described with its dimension (see
:mod:`ebbef2p.dimensions`).

Quick plots
-----------

.. code-block:: python

    from ebbef2p.plotting import plot_solution
    plot_solution(solution)
