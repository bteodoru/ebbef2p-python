Verification
============

Every case below is part of the automated test suite (``tests`` folder).

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
   * - Continuous beam on a two-parameter foundation with point load,
       concentrated moment and partial triangular load
     - Razaqpur & Shah (1991)
     - deflection 0.1 %; moment 52.18 vs 52.10 (0.15 %)
   * - Free beam on modified Vlasov foundation: central load, end loads,
       uniform load
     - Vallabhan & Das (1991): :math:`\gamma` = 0.507, 1.239, 0.395
     - 0.506, 1.240, 0.395
   * - Closed-form Vlasov parameters, homogeneous and linearly varying
       :math:`E`
     - numerical quadrature
     - :math:`10^{-10}`
   * - Global vertical equilibrium with variable :math:`k`, a jump in
       :math:`g`, supports, and soil beyond the ends
     - applied loads
     - :math:`10^{-8}`
   * - Beam shear force against :math:`\mathrm{d}M/\mathrm{d}x` with a
       shear layer
     - numerical differentiation
     - 0.2 %

The bending moment in the Razaqpur & Shah case is computed from element
equilibrium and does not change with the mesh; it is 0.15 % above the
published value. Moments obtained from the flexural stiffness alone
converge to the same 52.18 as the mesh is refined.
