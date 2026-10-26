"""Every numeric field of the model and the Vlasov classes declares its physical dimension."""

import dataclasses
import typing

import pytest

from ebbef2p.dimensions import (
    DIMENSIONLESS,
    FORCE,
    FORCE_PER_LENGTH,
    LENGTH,
    MOMENT,
    SECOND_MOMENT,
    STRESS,
)
from ebbef2p.model import (
    BeamSegment,
    DistributedLoad,
    FoundationSegment,
    PointLoad,
    Support,
)
from ebbef2p.vlasov import VlasovResult, VlasovSoil

# The table of HANDOFF stage 0, checked against docs/source/user_guide/conventions.rst:
# k and g are per unit length of beam, so k is F/L^2 and g is F.
EXPECTED = {
    BeamSegment: {"x0": LENGTH, "x1": LENGTH, "E": STRESS, "I": SECOND_MOMENT},
    FoundationSegment: {"x0": LENGTH, "x1": LENGTH, "k": STRESS, "g": FORCE},
    PointLoad: {"x": LENGTH, "force": FORCE, "moment": MOMENT},
    DistributedLoad: {"x0": LENGTH, "x1": LENGTH, "q0": FORCE_PER_LENGTH,
                      "q1": FORCE_PER_LENGTH},
    Support: {"x": LENGTH, "w": LENGTH, "theta": DIMENSIONLESS},
    VlasovSoil: {"E_top": STRESS, "nu": DIMENSIONLESS, "depth": LENGTH, "E_bottom": STRESS},
    VlasovResult: {"gamma": DIMENSIONLESS, "k": STRESS, "g": FORCE, "history": DIMENSIONLESS},
}


def is_numeric(hint) -> bool:
    return hint is float or float in typing.get_args(hint) or any(
        is_numeric(a) for a in typing.get_args(hint))


@pytest.mark.parametrize("cls", list(EXPECTED), ids=lambda c: c.__name__)
def test_every_numeric_field_has_dimension_and_doc(cls):
    hints = typing.get_type_hints(cls)
    numeric = {f.name for f in dataclasses.fields(cls) if is_numeric(hints[f.name])}
    assert numeric == set(EXPECTED[cls])
    for f in dataclasses.fields(cls):
        if f.name in numeric:
            assert f.metadata["doc"]


@pytest.mark.parametrize("cls", list(EXPECTED), ids=lambda c: c.__name__)
def test_dimension_table(cls):
    dims = {f.name: f.metadata.get("dimension") for f in dataclasses.fields(cls)}
    for name, dim in EXPECTED[cls].items():
        assert dims[name] == dim, name


def test_dimension_str():
    assert str(STRESS) == "F L^-2"
    assert str(MOMENT) == "F L"
    assert str(LENGTH) == "L"
    assert str(FORCE) == "F"
    assert str(DIMENSIONLESS) == "1"
