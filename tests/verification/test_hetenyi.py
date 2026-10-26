"""Long beam on a Winkler foundation (Hetenyi, 1946)."""

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
