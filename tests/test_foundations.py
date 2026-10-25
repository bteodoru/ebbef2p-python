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


# Razaqpur & Shah (1991), Int. J. Solids Struct. 27(4), 435-454, Example 2:
# continuous beam on a two-parameter foundation (EI = 2000, k = 64,
# k1 = g = 800), clamped at A, simply supported at B, C, D.

def razaqpur_example_2():
    L = 15.0
    m = Model().add_beam(0, L, E=2000, I=1)
    m.add_foundation(0, L, k=64, g=800)
    m.add_point_load(7.5, force=80)
    m.add_point_load(2, moment=100)
    m.add_distributed_load(10, L, 0, 18)
    m.add_support(0, w=0, theta=0)
    for x in (5, 10, L):
        m.add_support(x)
    return m


def test_razaqpur_exact_element_stiffness():
    """The exact 5 m element matrix printed in the paper is recovered by
    static condensation of a finely meshed 5 m segment."""
    from ebbef2p.mesh import build_mesh
    from ebbef2p.solver import assemble

    m = Model().add_beam(0, 5, E=2000, I=1).add_foundation(0, 5, k=64, g=800)
    mesh = build_mesh(m, 200)
    K = assemble(m, mesh)[0].toarray()
    n = K.shape[0]
    b = [0, 1, n - 2, n - 1]
    i = [d for d in range(n) if d not in b]
    S = K[np.ix_(b, b)] - K[np.ix_(b, i)] @ np.linalg.solve(K[np.ix_(i, i)], K[np.ix_(i, b)])

    paper = np.array([
        [498.05, 626.62, -341.54, 514.45],
        [626.62, 2131.05, -514.45, 659.46],
        [-341.54, -514.45, 498.05, -626.62],
        [514.45, 659.46, -626.62, 2131.05],
    ])
    assert S == pytest.approx(paper, abs=0.006)


def test_razaqpur_joint_rotations():
    """Exact joint rotations at B, C, D given in the paper (x 1e-3)."""
    s = solve(razaqpur_example_2(), n_elements=100)
    theta = [s.theta[s.mesh.node_index(x)] for x in (5, 10, 15)]
    assert theta == pytest.approx([4.249e-3, -6.096e-3, -7.157e-3], abs=6e-7)


def test_razaqpur_moment_diagram():
    """The paper gives the moment diagram only as a figure; at x = 2 the
    moment jumps by the applied 100 kNm, from about -48 to about +52."""
    s = solve(razaqpur_example_2(), n_elements=100)
    M = s.bending_moment()
    assert M.at(2, 'right') - M.at(2, 'left') == pytest.approx(100, rel=1e-9)
    assert M.at(2, 'right') == pytest.approx(52.18, abs=0.01)
    assert M.max() == pytest.approx(M.at(2, 'right'))
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
