Getting started
===============

ebbef2p-python is a Python implementation of the 1D Finite Element method for beams on elastic foundation. It allows you to do structural
analysis of beams on elastic foundation. It helps you to compute the forces and displacements in the structural elements.


Structure object
----------------

You start a model by instantiating a Structure object. All the models state, i.e. elements, materials and forces
are kept by this object.

.. autoclass:: ebbef2p.structure.Structure
