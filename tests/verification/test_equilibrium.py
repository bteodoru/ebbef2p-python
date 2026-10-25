"""Global vertical equilibrium: applied loads vs soil and support reactions."""

import pytest

from ebbef2p import Model, solve

L, E, I, q, P = 100.0, 432000.0, 2.25, 3.0, 20.0


def loaded_beam(**foundation):
    m = Model().add_beam(0, L, E, I)
    m.add_distributed_load(0, L, q).add_point_load(L/2, force=P)
    m.add_foundation(0, L, **foundation)
    return m


@pytest.mark.parametrize("soil_beyond_ends", [False, True])
def test_equilibrium_uniform_parameters(soil_beyond_ends):
    m = loaded_beam(k=6.4, g=700)
    m.soil_beyond_ends = soil_beyond_ends
    s = solve(m)
    assert s.reactive_pressure().resultant == pytest.approx(q*L + P, rel=1e-8)


def test_equilibrium_with_jump_in_g():
    m = Model().add_beam(0, L, E, I)
    m.add_distributed_load(0, L, q).add_point_load(L/2, force=P)
    m.add_foundation(0, L, k=(4, 9))
    m.add_foundation(0, 40, g=300).add_foundation(40, L, g=900)
    m.soil_beyond_ends = True
    s = solve(m)

    r = s.reactive_pressure()
    assert r.resultant == pytest.approx(q*L + P, rel=1e-8)
    assert any(x == pytest.approx(40) for x, _ in r.concentrated[1:-1])


def test_equilibrium_with_supports():
    m = loaded_beam(k=2.0, g=100)
    m.add_support(0).add_support(L)
    eq = solve(m).equilibrium()
    assert eq['supports'] > 0
    assert eq['residual'] == pytest.approx(0, abs=1e-7)
