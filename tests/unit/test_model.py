"""Model definition, validation and the inputs of the solvers."""

import numpy as np
import pytest

from ebbef2p import (
    ConvergenceError,
    EBBEF2pError,
    Model,
    ModelError,
    UnstableModelError,
    VlasovSoil,
    solve,
    solve_vlasov,
)

L = 100.0
SOIL = VlasovSoil(E_top=144.0, nu=0.2, depth=25.0)


def center_load(soil_beyond_ends=False):
    m = Model(soil_beyond_ends=soil_beyond_ends)
    return m.add_beam(0, 100, 432000.0, 2.25).add_point_load(50, force=20)


def test_singular_model_is_reported():
    m = Model().add_beam(0, 10, 1, 1).add_point_load(5, force=1)
    with pytest.raises(UnstableModelError):
        solve(m)


@pytest.mark.parametrize("build, message", [
    (lambda m: m.add_beam(0, 4, 1, 1).add_beam(5, 10, 1, 1), "gap"),
    (lambda m: m.add_beam(0, 10, 1, 1).add_point_load(12, force=1), "outside"),
])
def test_model_validation(build, message):
    m = Model()
    build(m)
    with pytest.raises(ModelError, match=message):
        solve(m)


def test_model_with_foundation_is_rejected():
    m = center_load().add_foundation(0, L, k=1)
    with pytest.raises(ModelError, match="foundation"):
        solve_vlasov(m, SOIL)


def test_soil_beyond_ends_is_taken_from_the_model():
    for flag in (True, False):
        r = solve_vlasov(center_load(flag), SOIL)
        assert r.solution.model.soil_beyond_ends is flag


@pytest.mark.parametrize("build", [
    lambda m: m.add_beam(0, 10, -1, 1),
    lambda m: m.add_beam(0, 10, 1, 1).add_foundation(4, 2, k=1),
    lambda m: m.add_beam(0, 10, 1, 1).add_foundation(0, 10, k=-1),
    lambda m: m.add_beam(0, 10, 1, 1).add_support(5, w=None, theta=None),
])
def test_invalid_input_is_rejected_when_added(build):
    with pytest.raises(ModelError):
        build(Model())


def test_conflicting_supports_are_rejected():
    m = Model().add_beam(0, 10, 1, 1).add_support(0).add_support(0, w=0.01).add_support(10)
    with pytest.raises(ModelError, match="conflicting"):
        solve(m)


@pytest.mark.parametrize("kwargs", [dict(nu=0.5), dict(depth=0.0), dict(E_bottom=-1.0)])
def test_invalid_vlasov_soil_is_rejected(kwargs):
    with pytest.raises(ModelError):
        VlasovSoil(**{**dict(E_top=144.0, nu=0.2, depth=25.0), **kwargs})


def test_gamma_iteration_reports_non_convergence():
    with pytest.raises(ConvergenceError, match="did not converge"):
        solve_vlasov(center_load(), SOIL, max_iter=1)


def test_exceptions_keep_their_builtin_bases():
    """Code written against the built-in exceptions keeps working."""
    assert issubclass(ModelError, ValueError)
    assert issubclass(UnstableModelError, np.linalg.LinAlgError)
    assert issubclass(ConvergenceError, RuntimeError)
    for exc in (ModelError, UnstableModelError, ConvergenceError):
        assert issubclass(exc, EBBEF2pError)
