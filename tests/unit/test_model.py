"""Model definition, validation and the inputs of the solvers."""

import numpy as np
import pytest

from ebbef2p import Model, VlasovSoil, solve, solve_vlasov

L = 100.0
SOIL = VlasovSoil(E_top=144.0, nu=0.2, depth=25.0)


def center_load(soil_beyond_ends=False):
    m = Model(soil_beyond_ends=soil_beyond_ends)
    return m.add_beam(0, 100, 432000.0, 2.25).add_point_load(50, force=20)


def test_singular_model_is_reported():
    m = Model().add_beam(0, 10, 1, 1).add_point_load(5, force=1)
    with pytest.raises(np.linalg.LinAlgError):
        solve(m)


@pytest.mark.parametrize("build, message", [
    (lambda m: m.add_beam(0, 4, 1, 1).add_beam(5, 10, 1, 1), "gap"),
    (lambda m: m.add_beam(0, 10, 1, 1).add_point_load(12, force=1), "outside"),
])
def test_model_validation(build, message):
    m = Model()
    build(m)
    with pytest.raises(ValueError, match=message):
        solve(m)


def test_model_with_foundation_is_rejected():
    m = center_load().add_foundation(0, L, k=1)
    with pytest.raises(ValueError, match="foundation"):
        solve_vlasov(m, SOIL)


def test_soil_beyond_ends_is_taken_from_the_model():
    for flag in (True, False):
        r = solve_vlasov(center_load(flag), SOIL)
        assert r.solution.model.soil_beyond_ends is flag
