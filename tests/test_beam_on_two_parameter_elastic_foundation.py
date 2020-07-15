import pytest
import numpy as np
from ebbef2p import Beam, Structure, NodalLoad, DistributedLoad, NodalSupport, ElasticFoundation


def test_razaqpur():
    """
    Razaqpur, A. G., & Shah, K. R. (1991). 
    Exact analysis of beams on two-parameter elastic foundations. 
    International Journal of Solids and Structures, 27(4), 435-454.
    """

    L = 15
    P = 80
    M = 100
    q = 18
    E = 24000
    h = 1
    w = 1
    k = 64
    t = 800
    
  
    w_max = -0.0254 #max displacement
    M_max = -52.10 #max moment
    tolerance = 1e-3 #set  a tolerance of 0.1%

    #initialize structure object
    s = Structure('test')

    #add beam
    s.add_beam(Beam(coord=[0, L],  E=E, h=h, w=w))

    #nodal loads
    s.add_load(NodalLoad(P, 7.5, 'fz'))
    s.add_load(NodalLoad(M, 2, 'my'))

    #distributed loads
    s.add_load(DistributedLoad(value=(0, q), position=(10, L)))

    #nodal support
    s.add_nodal_support(NodalSupport(position=0, uz=0, ur=0))
    s.add_nodal_support(NodalSupport(position=5, uz=0))
    s.add_nodal_support(NodalSupport(position=10, uz=0))
    s.add_nodal_support(NodalSupport(position=L, uz=0))

    #elastic foundation
    s.add_elastic_foundation(ElasticFoundation([k, k], [0, L], 'k'))    
    s.add_elastic_foundation(ElasticFoundation([t, t], [0, L], 't')) 
    s.discretize(100)
    s.add_elements()
    s.solve()
   
    assert min(s.get_displacements()['vertical_displacements']) == pytest.approx(w_max, rel=tolerance)
    assert min(s.get_bending_moments()['values']) == pytest.approx(M_max, rel=tolerance)