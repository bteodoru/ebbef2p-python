Finite element modelling
========================

Governing equation
------------------

For all two-parameter foundation models the reactive pressure is

.. math::
   :label: soil_pressure

   p(x) = k(x)\,w(x) - \frac{\mathrm{d}}{\mathrm{d}x}\left(g(x)\,\frac{\mathrm{d}w}{\mathrm{d}x}\right),

where :math:`k` is the first (Winkler) parameter and :math:`g` the second
(shear) parameter, :math:`g = 2t` in Vlasov's notation. For constant
:math:`g` this reduces to the familiar :math:`p = k\,w - g\,w''`. An
Euler-Bernoulli beam of bending stiffness :math:`EI` under a load
:math:`q(x)` then satisfies

.. math::
   :label: governing

   \frac{\mathrm{d}^2}{\mathrm{d}x^2}\left(EI\,\frac{\mathrm{d}^2w}{\mathrm{d}x^2}\right)
   - \frac{\mathrm{d}}{\mathrm{d}x}\left(g\,\frac{\mathrm{d}w}{\mathrm{d}x}\right)
   + k\,w = q.

Element formulation
-------------------

Each two-node element of length :math:`l` has the degrees of freedom
:math:`\mathbf{d}_e = [w_i, \theta_i, w_j, \theta_j]^T`, with
:math:`\theta = \mathrm{d}w/\mathrm{d}x`, interpolated by cubic Hermite
functions :math:`\mathbf{N}(\xi)`, :math:`\xi = (x - x_i)/l`. Galerkin's
method applied to :eq:`governing` gives

.. math::

   \left(\mathbf{K}_b + \mathbf{K}_k + \mathbf{K}_g\right)\mathbf{d}_e = \mathbf{f}_e,

   \mathbf{K}_b = \int_0^l EI\,\mathbf{N}''^T\mathbf{N}''\,\mathrm{d}x,\qquad
   \mathbf{K}_k = \int_0^l k\,\mathbf{N}^T\mathbf{N}\,\mathrm{d}x,\qquad
   \mathbf{K}_g = \int_0^l g\,\mathbf{N}'^T\mathbf{N}'\,\mathrm{d}x,\qquad
   \mathbf{f}_e = \int_0^l q\,\mathbf{N}^T\,\mathrm{d}x.

Within an element :math:`k`, :math:`g` and :math:`q` vary linearly between
their nodal values :math:`(k_i, k_j)`, :math:`(g_i, g_j)`,
:math:`(q_i, q_j)`, and the integrals are evaluated exactly:

.. math::

   \mathbf{K}_k = \begin{bmatrix}
   \frac{l(10k_i+3k_j)}{35} & \frac{l^2(15k_i+7k_j)}{420} & \frac{9l(k_i+k_j)}{140} & -\frac{l^2(7k_i+6k_j)}{420}\\
   & \frac{l^3(5k_i+3k_j)}{840} & \frac{l^2(6k_i+7k_j)}{420} & -\frac{l^3(k_i+k_j)}{280}\\
   & & \frac{l(3k_i+10k_j)}{35} & -\frac{l^2(7k_i+15k_j)}{420}\\
   \text{sym.} & & & \frac{l^3(3k_i+5k_j)}{840}
   \end{bmatrix}

.. math::

   \mathbf{K}_g = \begin{bmatrix}
   \frac{3(g_i+g_j)}{5l} & \frac{g_j}{10} & -\frac{3(g_i+g_j)}{5l} & \frac{g_i}{10}\\
   & \frac{l(3g_i+g_j)}{30} & -\frac{g_j}{10} & -\frac{l(g_i+g_j)}{60}\\
   & & \frac{3(g_i+g_j)}{5l} & -\frac{g_i}{10}\\
   \text{sym.} & & & \frac{l(g_i+3g_j)}{30}
   \end{bmatrix}

.. math::

   \mathbf{f}_e = \left[\frac{l(7q_i+3q_j)}{20},\;
   \frac{l^2(3q_i+2q_j)}{60},\;
   \frac{l(3q_i+7q_j)}{20},\;
   -\frac{l^2(2q_i+3q_j)}{60}\right]^T

:math:`\mathbf{K}_b` is the standard flexural stiffness matrix. The global
system is assembled in sparse form; prescribed displacements and rotations
are imposed by eliminating the corresponding degrees of freedom.

Soil beyond the beam ends
-------------------------

If the foundation continues indefinitely beyond a beam end, the soil there
carries no load and its surface deflection satisfies
:math:`k\,w - g\,w'' = 0`, decaying as

.. math::

   w(s) = w_\mathrm{end}\,e^{-\alpha s},\qquad \alpha = \sqrt{k/g},

with :math:`s` the distance from the end. The shear layer transmits to the
beam end a concentrated force :math:`g\,|w'(0^+)| = \sqrt{k\,g}\;w_\mathrm{end}`,
so the outside soil is represented exactly by a spring of stiffness
:math:`\sqrt{k\,g}` at each end node, with :math:`k` and :math:`g` taken at
that end. This option is set with ``Model(soil_beyond_ends=True)``.

Internal forces
---------------

Bending moments and shear forces are obtained from the equilibrium of each
element rather than by differentiating the interpolated displacements. The
nodal forces acting on element :math:`e` are

.. math::

   \mathbf{F}_e = \left(\mathbf{K}_b + \mathbf{K}_k + \mathbf{K}_g\right)\mathbf{d}_e - \mathbf{f}_e,

and, with :math:`M = -EI\,w''` and :math:`V = \mathrm{d}M/\mathrm{d}x`,

.. math::

   M_i = F_2,\quad M_j = -F_4,\qquad
   V_i = -F_1 - g_i\,\theta_i,\quad V_j = F_3 - g_j\,\theta_j.

The terms :math:`g\,\theta` remove the force carried by the shear layer, so
:math:`V` is the shear force in the beam itself. Computed this way, the end
moments do not depend on the mesh for problems such as point loads on
beams without foundation, and converge much faster than moments obtained
from :math:`\mathbf{K}_b\,\mathbf{d}_e` alone when a foundation is present.

Reactive pressure
-----------------

The distributed pressure :eq:`soil_pressure` is evaluated at both ends of
every element, with :math:`w''` from the Hermite interpolation:

.. math::

   p = k\,w - g\,w'' - g'\,w',\qquad g' = (g_j - g_i)/l.

Wherever the shear layer is cut, it transmits a concentrated force in
addition to :math:`p`. Integrating :math:`g\,w'\,v'` by parts in the weak
form gives these forces explicitly:

* at the start of the beam, :math:`-g\,w'`; at the end, :math:`+g\,w'`;
* at an interior node where :math:`g` jumps, the jump of :math:`g\,w'`
  (:math:`w'` itself is continuous);
* if the soil continues beyond the ends, the spring force
  :math:`\sqrt{k\,g}\;w_\mathrm{end}` at each end.

The total soil reaction is therefore

.. math::

   R = \int_0^L p\,\mathrm{d}x + \sum F_\mathrm{conc} = \int_0^L k\,w\,\mathrm{d}x
   + \sqrt{k\,g}\big|_0\,w(0) + \sqrt{k\,g}\big|_L\,w(L),

which :meth:`~ebbef2p.results.Solution.equilibrium` compares with the
applied loads and the support reactions. Omitting the concentrated forces
violates global equilibrium; in the uniform-load Vlasov example of the user
guide the two end forces carry about 10 % of the load.

References
----------

* Vallabhan, C. V. G., & Das, Y. C. (1991). Modified Vlasov model for beams
  on elastic foundations. *Journal of Geotechnical Engineering*, 117(6),
  956–966.
* Razaqpur, A. G., & Shah, K. R. (1991). Exact analysis of beams on
  two-parameter elastic foundations. *International Journal of Solids and
  Structures*, 27(4), 435–454.
* Hetényi, M. (1946). *Beams on Elastic Foundation*. University of
  Michigan Press.
