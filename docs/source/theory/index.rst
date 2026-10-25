Introduction
============

Beams on elastic foundations are used extensively by geotechnical,
pavement and railway engineers for the design and analysis of foundations.
Most programs for this problem implement only the simplest subgrade model,
Winkler's hypothesis, in which the soil is a bed of independent springs.
Such a model ignores the shear interaction between adjacent soil columns:
the surface displacement is discontinuous at the edges of a loaded area and
the soil outside the beam does not contribute.

Two-parameter models (Filonenko-Borodich, Pasternak, Vlasov) remove this
shortcoming by adding a second parameter that couples neighbouring points.
All of them lead to the same governing equation and differ only in how the
parameters are defined. The modified Vlasov model (Vallabhan & Das, 1991)
goes one step further and derives both parameters from the elastic
properties and the depth of the soil layer, so they need not be chosen
arbitrarily.

EBBEF2p was first written as a MATLAB toolbox and is available here as a
Python package. It uses two-node beam elements with cubic Hermite
interpolation; the element stiffness matrix and load vector follow from
Galerkin's method, with the foundation contributing two additional element
matrices. Because both foundation parameters may vary linearly within each
element, EBBEF2p accounts consistently for an inhomogeneous bearing soil.
Any static loading and support condition, including prescribed
displacements, can be analysed.
