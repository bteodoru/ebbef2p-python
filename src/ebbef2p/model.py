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
from itertools import pairwise
from typing import Self

from .dimensions import (
    DIMENSIONLESS,
    FORCE,
    FORCE_PER_LENGTH,
    LENGTH,
    MOMENT,
    SECOND_MOMENT,
    STRESS,
    quantity,
)
from .exceptions import ModelError

Pair = tuple[float, float]


_DOWN = "positive downwards"


@dataclass(frozen=True)
class BeamSegment:
    """Beam portion ``[x0, x1]`` with constant bending stiffness."""
    x0: float = field(metadata=quantity(LENGTH, "start of the segment"))
    x1: float = field(metadata=quantity(LENGTH, "end of the segment"))
    E: float = field(metadata=quantity(STRESS, "Young's modulus of the beam"))
    I: float = field(metadata=quantity(SECOND_MOMENT, "second moment of area of the section"))

    @property
    def EI(self) -> float:
        return self.E*self.I


@dataclass(frozen=True)
class FoundationSegment:
    """Two-parameter foundation acting on ``[x0, x1]``.

    ``k`` is the first (Winkler) parameter and ``g`` the second (shear)
    parameter, both given as ``(value at x0, value at x1)`` and varying
    linearly in between. ``g`` equals 2t in Vlasov's notation and G_p in
    Pasternak's. Overlapping segments add up.
    """
    x0: float = field(metadata=quantity(LENGTH, "start of the segment"))
    x1: float = field(metadata=quantity(LENGTH, "end of the segment"))
    k: Pair = field(default=(0.0, 0.0), metadata=quantity(
        STRESS, "first (Winkler) parameter per unit length of beam, at (x0, x1)"))
    g: Pair = field(default=(0.0, 0.0), metadata=quantity(
        FORCE, "second (shear) parameter per unit length of beam, at (x0, x1); "
        "g = 2t (Vlasov) = G_p (Pasternak)"))


@dataclass(frozen=True)
class PointLoad:
    """Concentrated force and/or moment at ``x``."""
    x: float = field(metadata=quantity(LENGTH, "position of the load"))
    force: float = field(default=0.0, metadata=quantity(FORCE, f"concentrated force, {_DOWN}"))
    moment: float = field(default=0.0, metadata=quantity(
        MOMENT, "concentrated moment, positive when it does positive work on a positive theta"))


@dataclass(frozen=True)
class DistributedLoad:
    """Load per unit length on ``[x0, x1]``, varying linearly from q0 to q1."""
    x0: float = field(metadata=quantity(LENGTH, "start of the load"))
    x1: float = field(metadata=quantity(LENGTH, "end of the load"))
    q0: float = field(metadata=quantity(FORCE_PER_LENGTH, f"load per unit length at x0, {_DOWN}"))
    q1: float = field(metadata=quantity(FORCE_PER_LENGTH, f"load per unit length at x1, {_DOWN}"))


@dataclass(frozen=True)
class Support:
    """Prescribed displacement and/or rotation at ``x``; ``None`` leaves it free."""
    x: float = field(metadata=quantity(LENGTH, "position of the support"))
    w: float | None = field(default=0.0, metadata=quantity(
        LENGTH, f"prescribed displacement, {_DOWN}; null leaves it free"))
    theta: float | None = field(default=None, metadata=quantity(
        DIMENSIONLESS, "prescribed rotation dw/dx; null leaves it free"))


def _pair(value: float | Pair) -> Pair:
    if isinstance(value, (int, float)):
        return (float(value), float(value))
    a, b = value
    return (float(a), float(b))


def _optional_float(value: float | None) -> float | None:
    return None if value is None else float(value)


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
    beams: list[BeamSegment] = field(default_factory=list)
    foundation: list[FoundationSegment] = field(default_factory=list)
    point_loads: list[PointLoad] = field(default_factory=list)
    distributed_loads: list[DistributedLoad] = field(default_factory=list)
    supports: list[Support] = field(default_factory=list)
    soil_beyond_ends: bool = field(default=False, metadata={
        "doc": "the foundation continues beyond both beam ends (end springs sqrt(k*g))"})

    # -- building -------------------------------------------------------
    def add_beam(self, x0: float, x1: float, E: float, I: float) -> Self:
        """Add a beam segment with constant E and I."""
        if x1 <= x0:
            raise ModelError("beam segment must have x1 > x0")
        if E <= 0 or I <= 0:
            raise ModelError("E and I must be positive")
        self.beams.append(BeamSegment(float(x0), float(x1), float(E), float(I)))
        return self

    def add_foundation(self, x0: float, x1: float, k: float | Pair = 0.0,
                       g: float | Pair = 0.0) -> Self:
        """Add a foundation segment; ``k`` and ``g`` are scalars or (start, end) pairs."""
        if x1 <= x0:
            raise ModelError("foundation segment must have x1 > x0")
        k, g = _pair(k), _pair(g)
        if min(k) < 0 or min(g) < 0:
            raise ModelError("foundation parameters must be non-negative")
        self.foundation.append(FoundationSegment(float(x0), float(x1), k, g))
        return self

    def add_point_load(self, x: float, force: float = 0.0, moment: float = 0.0) -> Self:
        """Add a concentrated force and/or moment at ``x``."""
        self.point_loads.append(PointLoad(float(x), float(force), float(moment)))
        return self

    def add_distributed_load(self, x0: float, x1: float, q0: float,
                             q1: float | None = None) -> Self:
        """Add a load varying from ``q0`` at ``x0`` to ``q1`` at ``x1`` (uniform if q1 is None)."""
        if x1 <= x0:
            raise ModelError("distributed load must have x1 > x0")
        q1 = q0 if q1 is None else q1
        self.distributed_loads.append(DistributedLoad(float(x0), float(x1), float(q0), float(q1)))
        return self

    def add_support(self, x: float, w: float | None = 0.0,
                    theta: float | None = None) -> Self:
        """Prescribe ``w`` and/or ``theta`` at ``x`` (``None`` = free).

        Defaults to a simple support (w = 0, rotation free). Use
        ``theta=0`` for a clamped end.
        """
        if w is None and theta is None:
            raise ModelError("a support must constrain w, theta or both")
        self.supports.append(Support(float(x), _optional_float(w), _optional_float(theta)))
        return self

    # -- queries ----------------------------------------------------------
    @property
    def extent(self) -> Pair:
        """(x_start, x_end) of the beam."""
        if not self.beams:
            raise ModelError("the model has no beam")
        return (min(b.x0 for b in self.beams), max(b.x1 for b in self.beams))

    def key_points(self) -> list[float]:
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

    def validate(self, tol: float = 1e-9) -> None:
        """Check that beam segments tile the beam and everything lies on it."""
        if not self.beams:
            raise ModelError("the model has no beam")
        segs = sorted(self.beams, key=lambda b: b.x0)
        for a, b in pairwise(segs):
            if abs(a.x1 - b.x0) > tol:
                raise ModelError(f"beam segments leave a gap or overlap at x = {a.x1}")
        x0, x1 = self.extent
        for x in self.key_points():
            if x < x0 - tol or x > x1 + tol:
                raise ModelError(f"x = {x} lies outside the beam [{x0}, {x1}]")
