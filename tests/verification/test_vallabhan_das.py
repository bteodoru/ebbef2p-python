"""Modified Vlasov model against Vallabhan & Das (1991)."""

import pytest

from ebbef2p import Model, VlasovSoil, solve_vlasov

# Vallabhan & Das (1991), J. Geotech. Eng. 117(6), 956-966:
# free beam, L = 100, H = L/4, Es = E/3000, nu = 0.2,
# soil continuing beyond the ends. Gamma reported to three decimals.
L, E, I = 100.0, 432000.0, 2.25
SOIL = VlasovSoil(E_top=E/3000, nu=0.2, depth=L/4)


def center_load(soil_beyond_ends=True):
    m = Model(soil_beyond_ends=soil_beyond_ends)
    return m.add_beam(0, L, E, I).add_point_load(L/2, force=20)


def end_loads(soil_beyond_ends=True):
    m = Model(soil_beyond_ends=soil_beyond_ends)
    return m.add_beam(0, L, E, I).add_point_load(0, force=20).add_point_load(L, force=20)


def uniform_load(soil_beyond_ends=True):
    m = Model(soil_beyond_ends=soil_beyond_ends)
    return m.add_beam(0, L, E, I).add_distributed_load(0, L, 3)


@pytest.mark.parametrize("model, gamma_ref", [
    (center_load, 0.507),
    (end_loads, 1.239),
    (uniform_load, 0.395),
])
def test_gamma_matches_reference(model, gamma_ref):
    r = solve_vlasov(model(), SOIL)
    assert r.gamma == pytest.approx(gamma_ref, abs=1.5e-3)
    assert r.solution.equilibrium()['residual'] == pytest.approx(0, abs=1e-6)


def test_soil_beyond_ends_matters_for_uniform_load():
    """Without the soil beyond the ends, a free beam under uniform load
    settles almost rigidly; with it, the ends settle markedly less."""
    ratio = {}
    for flag in (True, False):
        s = solve_vlasov(uniform_load(flag), SOIL).solution
        ratio[flag] = s.w[0]/s.w[len(s.w)//2]
    assert ratio[False] == pytest.approx(1.0, abs=0.01)
    assert ratio[True] < 0.8


def test_input_model_is_not_modified():
    m = center_load()
    solve_vlasov(m, SOIL)
    assert m.foundation == [] and m.soil_beyond_ends is True
