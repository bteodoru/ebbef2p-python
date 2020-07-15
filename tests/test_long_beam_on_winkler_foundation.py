import pytest
import numpy as np
from ebbef2p import Beam, Structure, NodalLoad, ElasticFoundation

def test_long_beam():

    L = 2
    P = 100
    E = 10
    h = 1
    w = 1
    I = w * h**3 / 12
    k = 10000
    
    characteristic_coefficient = (k/4/E/I)**0.25
    w_max = -P*characteristic_coefficient/2/k/w #analytical solution for max deflection
    tolerance = 1e-5 #set  a tolerance of 0.001%

    s = Structure('test')

    s.add_beam(Beam(coord=[0, L],  E=E, h=h, w=w))
    s.add_load(NodalLoad(P, L/2, 'fz'))
    s.add_elastic_foundation(ElasticFoundation([k, k], [0, L], 'k'))    
    s.add_elastic_foundation(ElasticFoundation([0, 0], [0, L], 't'))   
    s.discretize(100)
    s.add_elements()
    s.solve()
  
    #[np.where(s.nodes == 1)[0][0]]
    assert min(s.get_displacements()['vertical_displacements']) == pytest.approx(w_max, rel=tolerance)
