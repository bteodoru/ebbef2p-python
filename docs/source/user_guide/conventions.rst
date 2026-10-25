Conventions
===========

Units
    None are built in. With lengths in m and forces in kN, for example,
    ``E`` is in kPa, ``I`` in m\ :sup:`4`, distributed loads in kN/m,
    ``k`` in kN/m\ :sup:`2` (force per unit length of beam per unit
    displacement) and ``g`` in kN.

Coordinates and displacements
    ``x`` runs along the beam axis. The displacement ``w`` is positive
    downwards, so a settlement is positive, and
    :math:`\theta = \mathrm{d}w/\mathrm{d}x`.

Loads
    Forces and distributed loads are positive downwards. A concentrated
    moment is positive when it does positive work on a positive
    :math:`\theta`.

Internal forces
    The bending moment is :math:`M = -EI\,w''`, positive when the beam
    sags. The shear force is :math:`V = \mathrm{d}M/\mathrm{d}x` and is the
    force carried by the beam alone; the part carried by the foundation's
    shear layer, :math:`g\,w'`, is not included.

Reactions
    The soil pressure and the support reactions are positive when they act
    upwards on the beam, i.e. against positive loads, so that
    ``applied = soil + supports`` in
    :meth:`~ebbef2p.results.Solution.equilibrium`.

Foundation parameters
    ``k`` is the first (Winkler) parameter and ``g`` the second (shear)
    parameter. ``g`` equals :math:`2t` in Vlasov's notation and
    :math:`G_p` in Pasternak's. Both are per unit length of beam; for a
    beam of width :math:`b` on a soil characterised per unit area, multiply
    both by :math:`b`.
