"""Modified Vlasov model."""

import numpy as np
import pytest
from scipy.integrate import quad

from ebbef2p import Model, VlasovSoil, solve_vlasov

# Vallabhan & Das (1991), J. Geotech. Eng. 117(6), 956-966:
# free beam, L = 100, H = L/4, Es = E/3000, nu = 0.2,
# soil continuing beyond the ends. Gamma reported to three decimals.
L, E, I = 100.0, 432000.0, 2.25
SOIL = VlasovSoil(E_top=E/3000, nu=0.2, depth=L/4)


def center_load():
    return Model().add_beam(0, L, E, I).add_point_load(L/2, force=20)


def end_loads():
    return Model().add_beam(0, L, E, I).add_point_load(0, force=20).add_point_load(L, force=20)


def uniform_load():
    return Model().add_beam(0, L, E, I).add_distributed_load(0, L, 3)


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
        s = solve_vlasov(uniform_load(), SOIL, soil_beyond_ends=flag).solution
        ratio[flag] = s.w[0]/s.w[len(s.w)//2]
    assert ratio[False] == pytest.approx(1.0, abs=0.01)
    assert ratio[True] < 0.8


@pytest.mark.parametrize("gamma", [0.3, 1.0, 2.5])
@pytest.mark.parametrize("E_bottom", [None, 400.0])
def test_closed_form_parameters_match_quadrature(gamma, E_bottom):
    soil = VlasovSoil(E_top=150.0, nu=0.3, depth=12.0, E_bottom=E_bottom)
    E1 = soil.E_top
    E2 = E1 if E_bottom is None else E_bottom
    H, nu = soil.depth, soil.nu

    def E_z(z):
        return E1 + (E2 - E1)*z/H

    def phi(z):
        return np.sinh(gamma*(1 - z/H))/np.sinh(gamma)

    def dphi(z):
        return -gamma/H*np.cosh(gamma*(1 - z/H))/np.sinh(gamma)

    k_ref = (1 - nu)/((1 + nu)*(1 - 2*nu))*quad(lambda z: E_z(z)*dphi(z)**2, 0, H)[0]
    g_ref = 1/(2*(1 + nu))*quad(lambda z: E_z(z)*phi(z)**2, 0, H)[0]
    k, g = soil.parameters(gamma)
    assert k == pytest.approx(k_ref, rel=1e-10)
    assert g == pytest.approx(g_ref, rel=1e-10)


def test_model_with_foundation_is_rejected():
    m = center_load().add_foundation(0, L, k=1)
    with pytest.raises(ValueError, match="foundation"):
        solve_vlasov(m, SOIL)


def test_input_model_is_not_modified():
    m = center_load()
    solve_vlasov(m, SOIL)
    assert m.foundation == [] and m.soil_beyond_ends is False
