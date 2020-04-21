Finite Element Modelling
============================

All foundation models shown foregoing lead to the same differential equation.
Basically, all these models are equivalent and differ only in the definition of its
parameters [6]. 

 Governing Differential Equation
 -------------------------------

The various two-parameter elastic foundation models define the reactive pressure
of the foundation  :math:`p(x)`, as [6] 

.. math::
    :label: 2p_soil_pressure
    
    p(x) = kw(x) - t\frac{d^2w(x)}{dx^2}