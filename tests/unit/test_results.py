"""Post-processing: interpolation, consistency between results, mesh convergence."""

import numpy as np
import pytest

from ebbef2p import Model, solve

L, E, I, q, P = 100.0, 432000.0, 2.25, 3.0, 20.0


def loaded_beam(**foundation):
    m = Model().add_beam(0, L, E, I)
    m.add_distributed_load(0, L, q).add_point_load(L/2, force=P)
    m.add_foundation(0, L, **foundation)
    return m


def test_deflection_interpolates_between_nodes():
    m = Model().add_beam(0, 10, 1e4, 1).add_foundation(0, 10, k=100, g=50)
    m.add_distributed_load(0, 10, 5)
    s = solve(m, n_elements=10)
    assert s.deflection(s.x) == pytest.approx(s.w)
    mid = s.deflection([0.5, 3.3])
    assert np.all(np.isfinite(mid)) and mid.shape == (2,)


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


def test_winkler_pressure_is_k_times_w():
    s = solve(loaded_beam(k=6.4))
    r = s.reactive_pressure()
    assert r.values[0::2] == pytest.approx(6.4*s.w[:-1])
    assert all(f == pytest.approx(0) for _, f in r.concentrated)


def test_pressure_converges_with_mesh():
    m = loaded_beam(k=6.4, g=700)
    coarse, fine = solve(m, 50).reactive_pressure(), solve(m, 200).reactive_pressure()
    for x in (0.0, 20.0, 60.0, 100.0):
        assert coarse.at(x) == pytest.approx(fine.at(x), rel=1e-3)


def razaqpur_example_2():
    """Razaqpur & Shah (1991), example 2: a moment of 100 at x = 2."""
    m = Model().add_beam(0, 15, E=2000, I=1).add_foundation(0, 15, k=64, g=800)
    m.add_point_load(7.5, force=80).add_point_load(2, moment=100)
    m.add_distributed_load(10, 15, 0, 18)
    m.add_support(0, w=0, theta=0)
    for x in (5, 10, 15):
        m.add_support(x)
    return m


@pytest.fixture(scope="module")
def razaqpur():
    return solve(razaqpur_example_2(), n_elements=60)


def test_interpolate_matches_at_nodes(razaqpur):
    s = razaqpur
    for result in (s.bending_moment(), s.shear_force(), s.reactive_pressure()):
        for x in s.x[1:-1]:
            for side in ("left", "right"):
                assert result.interpolate(x, side) == result.at(x, side)


def test_interpolate_is_linear_inside_element(razaqpur):
    M = razaqpur.bending_moment()
    e = 7
    mid = 0.5*(M.x[2*e] + M.x[2*e + 1])
    value = M.interpolate(mid)
    assert isinstance(value, float)
    assert value == pytest.approx(0.5*(M.values[2*e] + M.values[2*e + 1]))
    many = M.interpolate(np.array([[mid, 1.0], [3.3, 14.9]]))
    assert many.shape == (2, 2) and many[0, 0] == pytest.approx(value)


def test_interpolate_keeps_moment_jump(razaqpur):
    M = razaqpur.bending_moment()
    jump = M.interpolate(2, "right") - M.interpolate(2, "left")
    assert jump == pytest.approx(M.at(2, "right") - M.at(2, "left"))
    assert abs(jump) == pytest.approx(100)


def test_interpolate_snaps_to_nodes(razaqpur):
    M = razaqpur.bending_moment()
    assert M.interpolate(2 + 1e-13, "left") == M.at(2, "left")
    assert M.interpolate(2 - 1e-13, "right") == M.at(2, "right")


def test_interpolate_outside_beam(razaqpur):
    M = razaqpur.bending_moment()
    for x in (-0.1, 15.1):
        with pytest.raises(ValueError, match="outside"):
            M.interpolate(x)
    for x in (0, 15):
        assert M.interpolate(x, "left") == M.interpolate(x, "right")
    with pytest.raises(ValueError, match="side"):
        M.interpolate(3, "middle")
    for x in (float("nan"), float("inf"), [1.0, float("nan")]):
        with pytest.raises(ValueError, match="finite"):
            M.interpolate(x)
