"""Vlasov foundation parameters against quadrature and their gamma -> 0 limits."""

import numpy as np
import pytest
from scipy.integrate import quad

from ebbef2p import Model, VlasovSoil, solve_vlasov


@pytest.mark.parametrize("gamma", [1e-8, 1e-4, 1e-2, 0.1, 0.29, 0.3, 1.0, 2.5])
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

    def integral(f):
        return quad(f, 0, H, epsabs=0, epsrel=1e-13)[0]

    k_ref = (1 - nu)/((1 + nu)*(1 - 2*nu))*integral(lambda z: E_z(z)*dphi(z)**2)
    g_ref = 1/(2*(1 + nu))*integral(lambda z: E_z(z)*phi(z)**2)
    k, g = soil.parameters(gamma)
    assert k == pytest.approx(k_ref, rel=1e-12)
    assert g == pytest.approx(g_ref, rel=1e-12)


def test_parameters_tend_to_linear_decay_as_gamma_vanishes():
    """For gamma -> 0, phi(z) = 1 - z/H and the parameters have closed
    limits: k = E_oed/H and g = E*H/(6(1+nu)) for a homogeneous layer."""
    soil = VlasovSoil(E_top=144.0, nu=0.2, depth=5.0)
    E_oed = 144.0*0.8/(1.2*0.6)
    for gamma in (0.0, 1e-12):
        k, g = soil.parameters(gamma)
        assert k == pytest.approx(E_oed/5.0, rel=1e-14)
        assert g == pytest.approx(144.0*5.0/(6*1.2), rel=1e-14)


def test_rigid_settlement_gives_vanishing_gamma():
    """A free beam under uniform load with no soil beyond the ends settles
    uniformly, so w' = 0, gamma -> 0 and g must stay finite and correct."""
    m = Model().add_beam(0, 10, E=30e6, I=0.02).add_distributed_load(0, 10, 10)
    r = solve_vlasov(m, VlasovSoil(E_top=144.0, nu=0.2, depth=5.0))
    assert r.gamma < 1e-6
    assert r.g == pytest.approx(100.0, rel=1e-10)
    assert r.solution.w == pytest.approx(10/r.k)
