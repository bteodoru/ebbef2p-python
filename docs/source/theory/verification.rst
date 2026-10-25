Verification
============

Every case below is part of the automated test suite, in
``tests/verification`` (except the last row, which is in ``tests/unit``).
Run them with ``uv run pytest tests/verification``.

.. list-table::
   :header-rows: 1
   :widths: 40 35 25

   * - Case
     - Reference
     - Agreement
   * - Element matrices :math:`\mathbf{K}_b`, :math:`\mathbf{K}_k`,
       :math:`\mathbf{K}_g`, :math:`\mathbf{f}_e` with linear :math:`k`,
       :math:`g`, :math:`q`
     - Gauss quadrature of the Galerkin integrals
     - round-off
   * - Simply supported beam, central point load and uniform load;
       cantilever; propped cantilever with prescribed settlement
     - closed-form beam theory
     - :math:`10^{-9}`
   * - Long beam on Winkler foundation, central point load
     - Hetényi (1946)
     - :math:`10^{-6}` (deflection), :math:`10^{-4}` (moment)
   * - Exact stiffness matrix of a 5 m beam element on a two-parameter
       foundation, recovered by static condensation of a fine mesh
     - Razaqpur & Shah (1991), exact element
     - all 16 entries to the printed precision
   * - Continuous beam on a two-parameter foundation with point load,
       concentrated moment and partial triangular load: joint rotations
     - Razaqpur & Shah (1991), exact solution
     - all three rotations to the four printed digits
   * - Free beam on modified Vlasov foundation: central load, end loads,
       uniform load
     - Vallabhan & Das (1991): :math:`\gamma` = 0.507, 1.239, 0.395
     - 0.506, 1.240, 0.395
   * - Vlasov parameters, homogeneous and linearly varying :math:`E`,
       :math:`10^{-8} \le \gamma \le 2.5` (closed forms above
       :math:`\gamma = 0.3`, Taylor series below)
     - numerical quadrature
     - :math:`10^{-12}`
   * - Vlasov parameters as :math:`\gamma \to 0`; free beam settling
       rigidly under uniform load
     - limits :math:`k = E_{oed}/H`, :math:`g = EH/(6(1+\nu))`;
       :math:`w = q/k`
     - :math:`10^{-14}`; :math:`10^{-10}`
   * - Global vertical equilibrium with variable :math:`k`, a jump in
       :math:`g`, supports, and soil beyond the ends
     - applied loads
     - :math:`10^{-8}`
   * - Beam shear force against :math:`\mathrm{d}M/\mathrm{d}x` with a
       shear layer
     - numerical differentiation
     - 0.2 %

Razaqpur & Shah give the internal forces of their Example 2 only as
diagrams. The moment computed here jumps from −47.82 to +52.18 under the
applied moment of 100 at x = 2, consistent with their figure; this value
does not change with the mesh because it comes from element equilibrium.
Note that their "vertical shear", :math:`-EI\,w''' + k_1 w'`, includes the
force carried by the shear layer, whereas
:meth:`~ebbef2p.results.Solution.shear_force` returns the shear force in
the beam alone, :math:`-EI\,w'''`.
