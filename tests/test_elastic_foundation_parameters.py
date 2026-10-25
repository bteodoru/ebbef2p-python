import pytest

from ebbef2p.structure import Structure
from ebbef2p.vlasov_foundation_parameters import VlasovFoundationParameters

"""
Girija Vallabhan, C. V., & Das, Y. C. (1991).
Modified Vlasov model for beams on elastic foundations.
Journal of geotechnical engineering, 117(6), 956-966.

Free beam, L = 100, H = L/4, Es = E/3000, nu = 0.2; the soil continues
beyond the beam ends. Reference gamma values are given to three decimals.
"""

L = 100
E = 432000
I = 2.25
nu = 0.2
P = 20
q = 3
EDEF = {'E1': E/3000, 'E2': E/3000}
H = L*0.25


def center_load(s):
    s.add_nodal_load(P, L/2, 'fz')


def end_loads(s):
    s.add_nodal_load(P, 0, 'fz')
    s.add_nodal_load(P, L, 'fz')


def uniform_load(s):
    s.add_distributed_load((q, q), (0, L))


@pytest.mark.parametrize("load, gamma_ref", [
    (center_load, 0.507),
    (end_loads, 1.239),
    (uniform_load, 0.395),
])
def test_modified_vlasov_gamma(load, gamma_ref):
    s = Structure('test')
    s.add_beam(coord=[0, L], E=E, I=I)
    load(s)
    s.add_nodes(50)

    p = s.compute_elastic_foundation_parameters(EDEF, nu, H)

    # gamma is reported to three decimals in the reference paper
    assert p[-1].gamma == pytest.approx(gamma_ref, abs=5e-4 + 1e-3)

    reference = VlasovFoundationParameters(EDEF, nu, H, gamma_ref)
    assert p[-1].k == pytest.approx(reference.k, rel=1e-3)
    assert p[-1].t == pytest.approx(reference.t, rel=1e-2)


def test_soil_beyond_ends_changes_uniform_load_response():
    """Without the soil beyond the ends a free beam under uniform load
    settles almost rigidly; with it, the ends settle markedly less."""
    results = {}
    for flag in (True, False):
        s = Structure('test')
        s.add_beam(coord=[0, L], E=E, I=I)
        uniform_load(s)
        s.add_nodes(50)
        s.compute_elastic_foundation_parameters(EDEF, nu, H, soil_beyond_ends=flag)
        w = s.get_displacements()['vertical_displacements']
        results[flag] = abs(w[0])/abs(w[len(w)//2])

    assert results[False] == pytest.approx(1.0, abs=0.01)
    assert results[True] < 0.8
