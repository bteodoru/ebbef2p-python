"""Razaqpur & Shah (1991): exact element stiffness and joint rotations."""

import numpy as np
import pytest

from ebbef2p import Model, solve


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
