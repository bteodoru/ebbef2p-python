"""Beams without a foundation, against closed-form solutions."""

import pytest

from ebbef2p import Model, solve

L, E, I = 2.0, 1.0, 1.0


def test_simply_supported_center_load():
    P = 100.0
    m = Model().add_beam(0, L, E, I).add_point_load(L/2, force=P)
    m.add_support(0).add_support(L)
    s = solve(m, n_elements=25)

    assert s.w.max() == pytest.approx(P*L**3/(48*E*I), rel=1e-9)
    assert s.bending_moment().max() == pytest.approx(P*L/4, rel=1e-9)
    V = s.shear_force()
    assert V.at(0, 'right') == pytest.approx(P/2, rel=1e-9)
    assert V.at(L, 'left') == pytest.approx(-P/2, rel=1e-9)
    assert [F for _, F, _ in s.support_reactions()] == pytest.approx([P/2, P/2])


def test_simply_supported_uniform_load():
    q = 10.0
    m = Model().add_beam(0, L, E, I).add_distributed_load(0, L, q)
    m.add_support(0).add_support(L)
    s = solve(m, n_elements=20)

    # Hermite elements are nodally exact for this case
    assert s.w.max() == pytest.approx(5*q*L**4/(384*E*I), rel=1e-9)
    assert s.bending_moment().max() == pytest.approx(q*L**2/8, rel=1e-9)
    assert s.shear_force().max() == pytest.approx(q*L/2, rel=1e-9)


def test_cantilever_uniform_load():
    q = 2.0
    m = Model().add_beam(0, L, E, I).add_distributed_load(0, L, q)
    m.add_support(0, w=0, theta=0)
    s = solve(m, n_elements=10)

    assert s.w[-1] == pytest.approx(q*L**4/(8*E*I), rel=1e-9)
    assert s.bending_moment().at(0, 'right') == pytest.approx(-q*L**2/2, rel=1e-9)
    assert s.shear_force().at(0, 'right') == pytest.approx(q*L, rel=1e-9)


def test_prescribed_settlement():
    """Propped cantilever with a settlement d at the prop: M(0) = 3EI d/L^2."""
    d = 0.01
    m = Model().add_beam(0, L, E, I)
    m.add_support(0, w=0, theta=0).add_support(L, w=d)
    s = solve(m, n_elements=10)
    assert abs(s.bending_moment().at(0, 'right')) == pytest.approx(3*E*I*d/L**2, rel=1e-9)


def test_stepped_beam_continuity():
    """Two segments with different EI: displacements stay continuous and the
    moment is continuous across the change of section."""
    m = Model().add_beam(0, 1, E, I).add_beam(1, 2, 4*E, I)
    m.add_point_load(1.5, force=10).add_support(0).add_support(2)
    s = solve(m, n_elements=40)
    M = s.bending_moment()
    assert M.at(1, 'left') == pytest.approx(M.at(1, 'right'), rel=1e-9)
    assert s.equilibrium()['residual'] == pytest.approx(0, abs=1e-8)
