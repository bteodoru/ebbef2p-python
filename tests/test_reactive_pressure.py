import numpy as np
import pytest

from ebbef2p.structure import Structure

L = 100
E = 432000
I = 2.25
q = 3
P = 20


def free_beam(k, t, n=100, soil_beyond_ends=False):
    s = Structure('test')
    s.add_beam(coord=[0, L], E=E, I=I)
    s.add_distributed_load((q, q), (0, L))
    s.add_nodal_load(P, L/2, 'fz')
    s.add_elastic_foundation(k[0], k[1], 'k')
    s.add_elastic_foundation(t[0], t[1], 't')
    s.soil_beyond_ends = soil_beyond_ends
    s.add_nodes(n)
    s.add_elements(s.nodes)
    s.solve(s.build_global_matrix(), s.build_load_vector(), s.get_boudary_conditions())
    return s


@pytest.mark.parametrize("soil_beyond_ends", [False, True])
def test_equilibrium_uniform_parameters(soil_beyond_ends):
    s = free_beam(([6.4, 6.4], [0, L]), ([700.0, 700.0], [0, L]),
                  soil_beyond_ends=soil_beyond_ends)
    r = s.get_reactive_pressure()
    assert r['resultant'] == pytest.approx(q*L + P, rel=1e-6)


def test_equilibrium_variable_parameters_with_jump_in_t():
    """k varies linearly; t jumps at x = 40 (two foundation zones)."""
    s = Structure('test')
    s.add_beam(coord=[0, L], E=E, I=I)
    s.add_distributed_load((q, q), (0, L))
    s.add_nodal_load(P, L/2, 'fz')
    s.add_elastic_foundation([4.0, 9.0], [0, L], 'k')
    s.add_elastic_foundation([300.0, 300.0], [0, 40], 't')
    s.add_elastic_foundation([900.0, 900.0], [40, L], 't')
    s.soil_beyond_ends = True
    s.add_nodes(100)
    s.add_elements(s.nodes)
    s.solve(s.build_global_matrix(), s.build_load_vector(), s.get_boudary_conditions())

    r = s.get_reactive_pressure()
    assert r['resultant'] == pytest.approx(q*L + P, rel=1e-6)
    # a concentrated force must appear where t jumps
    assert any(x == pytest.approx(40) for x, _ in r['concentrated'][1:-1])


def test_winkler_pressure_is_k_times_w():
    s = free_beam(([6.4, 6.4], [0, L]), ([0.0, 0.0], [0, L]))
    r = s.get_reactive_pressure()
    u = np.asarray(s.u).flatten()
    w_nodes = u[0::2]
    # values at the start of each element vs k*w at that node
    assert np.allclose(r['values'][0::2], 6.4*w_nodes[:-1])
    assert all(f == pytest.approx(0) for _, f in r['concentrated'])


def test_pressure_converges_with_mesh():
    coarse = free_beam(([6.4, 6.4], [0, L]), ([700.0, 700.0], [0, L]), n=50)
    fine = free_beam(([6.4, 6.4], [0, L]), ([700.0, 700.0], [0, L]), n=200)
    pc = dict(zip(coarse.get_reactive_pressure()['coords'], coarse.get_reactive_pressure()['values']))
    pf = dict(zip(fine.get_reactive_pressure()['coords'], fine.get_reactive_pressure()['values']))
    for x in (0.0, 20.0, 60.0, 100.0):
        assert pc[x] == pytest.approx(pf[x], rel=1e-3)
