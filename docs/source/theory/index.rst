Introduction
============

The concept of beams on elastic foundations it is extensively used by geotechnical,
pavement and railroad engineers for foundation design and analysis.

Currently, the analysis of beams on elastic foundation is performed by using
special computer programs based on numerical methods, such as Finite Difference
Method (FDM) and Finite Element Method (FEM). However, these programs are
limited in their application, most of them being developed only for a very simple
subgrade model, Winkler's Hypothesis. They cannot be used for other soil models
such as Two-Parameter, Elastic Half-Space or Elastic Layer and others.

This paper describes a finite element computer program, as a toolbox to MATLAB,
developed to analyse the interaction between a beam and its two-parameter elastic
foundation. By considering a linear variation of both foundation parameter,
EBBEF2p can account in a consistent way for the bearing soil inhomogeneity. It
can be used for any practical static loading and support condition including
prescribed displacement.

The numerical model uses a cubic Hermitian polynomial to interpolate nodal
values of the displacements field for a two-node beam elements. The elemental
stiffness matrix and load vector are obtained by using Galerkin’s Residual Method 
with adding the contribution of the foundation as element foundation stiffness
matrices to the regular flexure beam element. 