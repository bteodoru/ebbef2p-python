"""Beams on Winkler and two-parameter foundations."""

import numpy as np
import pytest

from ebbef2p import Model, solve


def test_long_beam_on_winkler_foundation():
    """Hetenyi: long beam, central point load, w_max = P*lambda/(2k)."""
    L, P, E, I, k = 2.0, 100.0, 1.0, 1.0, 10000.0
    lam = (k/(4*E*I))**0.25
    m = Model().add_beam(0, L, E, I).add_foundation(0, L, k=k)
    m.add_point_load(L/2, force=P)
    s = solve(m, n_elements=100)

    assert s.w.max() == pytest.approx(P*lam/(2*k), rel=1e-6)
    assert s.bending_moment().max() == pytest.approx(P/(4*lam), rel=1e-4)


def test_razaqpur_two_parameter():
    """Razaqpur & Shah (1991), Int. J. Solids Struct. 27(4), 435-454.

    Continuous beam on a two-parameter foundation with a point load, a
    concentrated moment and a partial uniform load.
    """
    L = 15.0
    m = Model().add_beam(0, L, E=2000, I=1)
    m.add_foundation(0, L, k=64, g=800)
    m.add_point_load(7.5, force=80)
    m.add_point_load(2, moment=100)
    m.add_distributed_load(10, L, 0, 18)
    m.add_support(0, w=0, theta=0)
    for x in (5, 10, L):
        m.add_support(x)
    s = solve(m, n_elements=100)

    assert s.w.max() == pytest.approx(0.0254, rel=1e-3)
    # The equilibrium-based moment is mesh-independent here and gives
    # 52.18, 0.15 % above the published 52.10.
    assert s.bending_moment().max() == pytest.approx(52.10, rel=2e-3)
    assert s.equilibrium()['residual'] == pytest.approx(0, abs=1e-8)


def test_beam_shear_is_derivative_of_moment():
    """With a shear layer, the beam shear (not the total transverse force)
    must satisfy V = dM/dx."""
    m = Model().add_beam(0, 20, 30e6, 0.01)
    m.add_foundation(0, 20, k=(4000, 9000), g=(2000, 5000))
    m.add_point_load(8, force=300).add_distributed_load(10, 20, 20, 5)
    s = solve(m, n_elements=400)

    M, V = s.bending_moment(), s.shear_force()
    x, Mi = M.x[0::2], M.values[0::2]
    dMdx = np.gradient(Mi, x)
    for i in (50, 150, 300, 350):
        assert V.values[2*i] == pytest.approx(dMdx[i], rel=2e-3)


def test_results_converge_with_mesh():
    m = Model().add_beam(0, 20, 30e6, 0.01)
    m.add_foundation(0, 20, k=(4000, 9000), g=(2000, 5000))
    m.add_point_load(8, force=300)
    coarse, fine = solve(m, 20), solve(m, 320)
    assert coarse.deflection([0, 8, 20]) == pytest.approx(fine.deflection([0, 8, 20]), rel=1e-4)
    assert coarse.bending_moment().at(8) == pytest.approx(fine.bending_moment().at(8), rel=1e-3)


def test_deflection_interpolates_between_nodes():
    m = Model().add_beam(0, 10, 1e4, 1).add_foundation(0, 10, k=100, g=50)
    m.add_distributed_load(0, 10, 5)
    s = solve(m, n_elements=10)
    assert s.deflection(s.x) == pytest.approx(s.w)
    mid = s.deflection([0.5, 3.3])
    assert np.all(np.isfinite(mid)) and mid.shape == (2,)


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
