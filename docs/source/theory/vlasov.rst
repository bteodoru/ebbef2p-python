Modified Vlasov model
=====================

Vlasov and Leont'ev treated the soil as an elastic layer of depth
:math:`H` on a rigid base and assumed the vertical displacement in the
layer to be :math:`w(x)\,\phi(z)`, with :math:`\phi(0) = 1` and
:math:`\phi(H) = 0`. Vallabhan & Das (1991) chose

.. math::

   \phi(z) = \frac{\sinh\left(\gamma\,(1 - z/H)\right)}{\sinh\gamma},

where the non-dimensional parameter :math:`\gamma` controls how fast the
displacement decays with depth.

Foundation parameters
---------------------

For a layer with Poisson's ratio :math:`\nu` and a Young's modulus varying
linearly from :math:`E_1` at the surface to :math:`E_2` at the base,

.. math::

   k = \frac{1-\nu}{(1+\nu)(1-2\nu)}\int_0^H E(z)\,\phi'(z)^2\,\mathrm{d}z,
   \qquad
   g = 2t = \frac{1}{2(1+\nu)}\int_0^H E(z)\,\phi(z)^2\,\mathrm{d}z,

which evaluate in closed form to

.. math::

   k = \frac{1-\nu}{(1+\nu)(1-2\nu)H}\,
   \frac{E_1\left(2\gamma\sinh 2\gamma + 4\gamma^2\right)
   + (E_2-E_1)\left(\cosh 2\gamma - 1 + 2\gamma^2\right)}{8\sinh^2\gamma},

.. math::

   g = \frac{H}{(1+\nu)\gamma^2}\,
   \frac{E_1\left(2\gamma\sinh 2\gamma - 4\gamma^2\right)
   + (E_2-E_1)\left(\cosh 2\gamma - 1 - 2\gamma^2\right)}{16\sinh^2\gamma}.

Both expressions are checked against numerical quadrature in the test
suite.

Iteration for gamma
-------------------

Minimising the total potential energy of the beam–soil system with respect
to :math:`\phi` gives

.. math::

   \left(\frac{\gamma}{H}\right)^2 = \frac{1-2\nu}{2(1-\nu)}\,
   \frac{\displaystyle\int_{-\infty}^{\infty}\left(\frac{\mathrm{d}w}{\mathrm{d}x}\right)^2\mathrm{d}x}
        {\displaystyle\int_{-\infty}^{\infty}w^2\,\mathrm{d}x}.

The integrals run over the whole soil surface. Inside the beam they are
evaluated by Gauss quadrature of the finite element interpolation. If the
soil continues beyond the ends, :math:`w = w_\mathrm{end}\,e^{-\alpha s}`
there, with :math:`\alpha = \sqrt{k/g}`, and the outside contributions are

.. math::

   \int w^2\,\mathrm{d}s = \frac{w_\mathrm{end}^2}{2\alpha},\qquad
   \int w'^2\,\mathrm{d}s = \frac{\alpha\,w_\mathrm{end}^2}{2}

at each end.

:func:`~ebbef2p.vlasov.solve_vlasov` starts from :math:`\gamma = 1`,
computes :math:`k` and :math:`g`, solves the beam, updates :math:`\gamma`
from the deflected shape, and repeats until :math:`\gamma` changes by less
than the tolerance. Convergence usually takes 5 to 10 iterations.

.. important::

   The assumption about the soil beyond the ends must be the same in the
   stiffness matrix (the end springs :math:`\sqrt{k\,g}`) and in the update
   of :math:`\gamma`. Using the outside terms in the :math:`\gamma` update
   while omitting the end springs gives a beam that settles almost rigidly
   under uniform load and a :math:`\gamma` about 26 % too high for the
   uniform-load case of Vallabhan & Das. Both are controlled together by
   ``Model(soil_beyond_ends=...)``; :func:`~ebbef2p.vlasov.solve_vlasov`
   reads the flag from the model and has no separate option for it.

Small values of :math:`\gamma`
-------------------------------

When the beam settles almost uniformly, :math:`w' \to 0` and
:math:`\gamma \to 0`; the decay function tends to :math:`\phi(z) = 1 - z/H`
and, for a homogeneous layer, :math:`k \to E(1-\nu)/((1+\nu)(1-2\nu)H)`
and :math:`g \to E H/(6(1+\nu))`. The closed forms above cancel badly in
this limit (the relative error of :math:`g` grows like
:math:`10^{-16}/\gamma^4`), so for :math:`\gamma < 0.3` their Taylor series
in :math:`\gamma^2`, up to :math:`\gamma^{14}`, is used instead. Both agree
to :math:`10^{-15}` at the switch.
