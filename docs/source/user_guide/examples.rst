Examples
========

The scripts below are in the ``examples`` folder of the repository and are
run by the test suite, so they always match the current version.

Long beam on a Winkler foundation
---------------------------------

A central point load on a long beam, compared with Hetényi's solution.

.. literalinclude:: ../../../examples/winkler_long_beam.py
   :language: python

Footing crossing a weak zone
----------------------------

Both foundation parameters vary along the footing: ``k`` drops over a soft
lens and ``g`` changes abruptly. The soil continues beyond the ends. The
concentrated soil forces appear at the ends and where ``g`` jumps.

.. literalinclude:: ../../../examples/variable_foundation.py
   :language: python

Modified Vlasov foundation
--------------------------

The uniform-load case of Vallabhan & Das (1991). The foundation parameters
follow from the elastic layer and the iterated decay parameter
:math:`\gamma`.

.. literalinclude:: ../../../examples/modified_vlasov.py
   :language: python
