"""Physical dimensions of the model fields.

The package has no built-in units; every numeric field of the model
dataclasses declares instead its dimension as exponents of force (F) and
length (L) in ``field(metadata={"dimension": ..., "doc": ...})``. With
lengths in m and forces in kN, for example, ``STRESS`` is kPa.

``k`` is a force per unit length of beam per unit displacement, so it has
the dimension of a stress (F L^-2) although it is not one; ``g`` is a force.
"""

from typing import Any, NamedTuple


class Dimension(NamedTuple):
    """Exponents of force and length, e.g. ``Dimension(1, -2)`` for F L^-2."""
    force: int
    length: int

    def __str__(self) -> str:
        parts = [f"{symbol}^{power}" if power != 1 else symbol
                 for symbol, power in (("F", self.force), ("L", self.length)) if power]
        return " ".join(parts) or "1"


LENGTH = Dimension(0, 1)
FORCE = Dimension(1, 0)
MOMENT = Dimension(1, 1)
STRESS = Dimension(1, -2)
FORCE_PER_LENGTH = Dimension(1, -1)
SECOND_MOMENT = Dimension(0, 4)
DIMENSIONLESS = Dimension(0, 0)


def quantity(dimension: Dimension, doc: str) -> dict[str, Any]:
    """Field metadata: ``field(metadata=quantity(LENGTH, "start of the segment"))``."""
    return {"dimension": dimension, "doc": doc}
