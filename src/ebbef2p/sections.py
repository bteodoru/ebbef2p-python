"""Section properties of the beam."""

from .exceptions import ModelError


def rectangular_section(b: float, h: float) -> float:
    """Second moment of area ``I = b*h**3/12`` of a rectangular section.

    Parameters
    ----------
    b : float
        Width of the beam in contact with the soil.
    h : float
        Height of the section.

    Notes
    -----
    For a section that is not rectangular (a T or inverted-T footing, a
    rail), enter the equivalent height ``h_eq = (12*I/b)**(1/3)`` computed
    with the contact width ``b``, not a width chosen only to reproduce I.

    The same ``b`` converts foundation parameters given per unit area into
    the per-unit-length values of the model: ``k = k_s*b`` (and
    ``g = g_s*b``). A width that matches I but not the contact with the soil
    would give a wrong ``k``.
    """
    if b <= 0 or h <= 0:
        raise ModelError("b and h must be positive")
    return b*h**3/12
