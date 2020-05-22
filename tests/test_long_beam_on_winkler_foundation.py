import pytest
import numpy as np
from ebbef2p.structure import Structure


def test_long_beam():

    L = 2
    P = 100
    E = 1
    I = 1
    k = 10000
    
    characteristic_coefficient = (k/4/E/I)**0.25
    w_max = -P*characteristic_coefficient/2/k #analytical solution for max deflection
    tolerance = 1e-6 #set  a tolerance of 0.0001%

    s = Structure('test')

    s.add_beam(coord=[0, L],  E=E, I=I)
    s.add_nodal_load(P, 1, 'fz')
    s.add_elastic_foundation([k, k], [0, 2], 'k')
    s.add_elastic_foundation([0, 0], [0, 2], 't')
    s.add_nodes(100)
    s.add_elements(s.nodes)
    s.solve(s.build_global_matrix(), s.build_load_vector(), s.get_boudary_conditions())
    #[np.where(s.nodes == 1)[0][0]]
    assert min(s.get_displacements()['vertical_displacements']) == pytest.approx(w_max, rel=tolerance)
