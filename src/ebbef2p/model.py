"""Problem definition: beam segments, foundation, loads and supports.

Conventions
-----------
* ``x`` runs along the beam axis; ``w`` is the transverse displacement,
  positive downwards (i.e. a settlement is positive).
* ``theta = dw/dx``.
* Forces and distributed loads are positive downwards. A concentrated moment
  is positive when it does positive work on a positive ``theta``.
* No system of units is built in; use any consistent set.
"""

from dataclasses import dataclass, field
from typing import List, Optional, Tuple


@dataclass(frozen=True)
class BeamSegment:
    """Beam portion ``[x0, x1]`` with constant bending stiffness."""
    x0: float
    x1: float
    E: float
    I: float

    @property
    def EI(self):
        return self.E*self.I


@dataclass(frozen=True)
class FoundationSegment:
    """Two-parameter foundation acting on ``[x0, x1]``.

    ``k`` is the first (Winkler) parameter and ``g`` the second (shear)
    parameter, both given as ``(value at x0, value at x1)`` and varying
    linearly in between. ``g`` equals 2t in Vlasov's notation and G_p in
    Pasternak's. Overlapping segments add up.
    """
    x0: float
    x1: float
    k: Tuple[float, float] = (0.0, 0.0)
    g: Tuple[float, float] = (0.0, 0.0)


@dataclass(frozen=True)
class PointLoad:
    x: float
    force: float = 0.0
    moment: float = 0.0


@dataclass(frozen=True)
class DistributedLoad:
    """Load per unit length on ``[x0, x1]``, varying linearly from q0 to q1."""
    x0: float
    x1: float
    q0: float
    q1: float


@dataclass(frozen=True)
class Support:
    """Prescribed displacement and/or rotation at ``x``; ``None`` leaves it free."""
    x: float
    w: Optional[float] = 0.0
    theta: Optional[float] = None


def _pair(value):
    if isinstance(value, (int, float)):
        return (float(value), float(value))
    a, b = value
    return (float(a), float(b))


@dataclass
class Model:
    """A beam on a two-parameter elastic foundation.

    Parameters
    ----------
    soil_beyond_ends : bool
        If True, the foundation continues indefinitely beyond both beam ends.
        Its effect is represented by an equivalent spring sqrt(k*g) at each
        end, with k and g taken at that end (Vallabhan & Das, 1991).

    Examples
    --------
    >>> m = Model()
    >>> m.add_beam(0, 10, E=30e6, I=0.01)
    >>> m.add_foundation(0, 10, k=5000, g=2000)
    >>> m.add_point_load(5, force=100)
    """
    beams: List[BeamSegment] = field(default_factory=list)
    foundation: List[FoundationSegment] = field(default_factory=list)
    point_loads: List[PointLoad] = field(default_factory=list)
    distributed_loads: List[DistributedLoad] = field(default_factory=list)
    supports: List[Support] = field(default_factory=list)
    soil_beyond_ends: bool = False

    # -- building -------------------------------------------------------
    def add_beam(self, x0, x1, E, I):
        """Add a beam segment with constant E and I."""
        if x1 <= x0:
            raise ValueError("beam segment must have x1 > x0")
        if E <= 0 or I <= 0:
            raise ValueError("E and I must be positive")
        self.beams.append(BeamSegment(float(x0), float(x1), float(E), float(I)))
        return self

    def add_foundation(self, x0, x1, k=0.0, g=0.0):
        """Add a foundation segment; ``k`` and ``g`` are scalars or (start, end) pairs."""
        if x1 <= x0:
            raise ValueError("foundation segment must have x1 > x0")
        k, g = _pair(k), _pair(g)
        if min(k) < 0 or min(g) < 0:
            raise ValueError("foundation parameters must be non-negative")
        self.foundation.append(FoundationSegment(float(x0), float(x1), k, g))
        return self

    def add_point_load(self, x, force=0.0, moment=0.0):
        """Add a concentrated force and/or moment at ``x``."""
        self.point_loads.append(PointLoad(float(x), float(force), float(moment)))
        return self

    def add_distributed_load(self, x0, x1, q0, q1=None):
        """Add a distributed load from ``q0`` at ``x0`` to ``q1`` at ``x1`` (uniform if q1 is None)."""
        if x1 <= x0:
            raise ValueError("distributed load must have x1 > x0")
        q1 = q0 if q1 is None else q1
        self.distributed_loads.append(DistributedLoad(float(x0), float(x1), float(q0), float(q1)))
        return self

    def add_support(self, x, w=0.0, theta=None):
        """Prescribe ``w`` and/or ``theta`` at ``x`` (``None`` = free).

        Defaults to a simple support (w = 0, rotation free). Use
        ``theta=0`` for a clamped end.
        """
        if w is None and theta is None:
            raise ValueError("a support must constrain w, theta or both")
        self.supports.append(Support(float(x), w, theta))
        return self

    # -- queries ----------------------------------------------------------
    @property
    def extent(self):
        """(x_start, x_end) of the beam."""
        if not self.beams:
            raise ValueError("the model has no beam")
        return (min(b.x0 for b in self.beams), max(b.x1 for b in self.beams))

    def key_points(self):
        """Coordinates that must coincide with mesh nodes."""
        pts = set()
        for b in self.beams:
            pts.update((b.x0, b.x1))
        for f in self.foundation:
            pts.update((f.x0, f.x1))
        for d in self.distributed_loads:
            pts.update((d.x0, d.x1))
        pts.update(p.x for p in self.point_loads)
        pts.update(s.x for s in self.supports)
        return sorted(pts)

    def validate(self, tol=1e-9):
        """Check that beam segments tile the beam and everything lies on it."""
        if not self.beams:
            raise ValueError("the model has no beam")
        segs = sorted(self.beams, key=lambda b: b.x0)
        for a, b in zip(segs, segs[1:]):
            if abs(a.x1 - b.x0) > tol:
                raise ValueError(f"beam segments leave a gap or overlap at x = {a.x1}")
        x0, x1 = self.extent
        for x in self.key_points():
            if x < x0 - tol or x > x1 + tol:
                raise ValueError(f"x = {x} lies outside the beam [{x0}, {x1}]")
