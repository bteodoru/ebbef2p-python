import pytest
import numpy as np
from ebbef2p.structure import Structure


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
    E = 2000
    I = 1
    k = 64
    t = 800
    
  
    w_max = -0.0254 #max displacement
    M_max = 52.10 #max moment
    tolerance = 1e-3 #set  a tolerance of 0.1%

    #initialize structure object
    s = Structure('test')

    #add beam
    s.add_beam(coord=[0, L],  E=E, I=I)

    #nodal loads
    s.add_nodal_load(P, 7.5, 'fz')
    s.add_nodal_load(M, 2, 'my')

    #distributed loads
    s.add_distributed_load((0, q), (10, L))

    #nodal support
    s.add_nodal_support({'uz': 0, 'ur': 0}, 0)
    s.add_nodal_support({'uz': 0, 'ur': "NaN"}, 5)
    s.add_nodal_support({'uz': 0, 'ur': "NaN"}, 10)
    s.add_nodal_support({'uz': 0, 'ur': "NaN"}, L)

    #elastic foundation
    s.add_elastic_foundation([k, k], [0, L], 'k')
    s.add_elastic_foundation([t, t], [0, L], 't')
    s.add_nodes(100)
    s.add_elements(s.nodes)
    s.solve(s.build_global_matrix(), s.build_load_vector(), s.get_boudary_conditions())
   
    assert min(s.get_displacements()['vertical_displacements']) == pytest.approx(w_max, rel=tolerance)
    assert max(s.get_bending_moments()['values']) == pytest.approx(M_max, rel=tolerance)