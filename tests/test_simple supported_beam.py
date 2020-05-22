import pytest
import numpy as np
from ebbef2p.structure import Structure

L = 2
E = 1
I = 1

def test_center_load():

    P = 100
       
    M_max = P * L / 4  # maximum moment
    S_max = P/2 # max shearing force
    w_max = -P * L ** 3 / (48 * E * I)  # max displacement
    tolerance = 1e-6 #set  a tolerance of 0.0001%

    s = Structure('test')

    s.add_beam(coord=[0, L],  E=E, I=I)
    s.add_nodal_load(P, L/2, 'fz')
    s.add_nodal_support({'uz': 0, 'ur': "NaN"}, 0)
    s.add_nodal_support({'uz': 0, 'ur': "NaN"}, L)
    s.add_nodes(25)
    s.add_elements(s.nodes)
    s.solve(s.build_global_matrix(), s.build_load_vector(), s.get_boudary_conditions())
  
    assert min(s.get_displacements()['vertical_displacements']) == pytest.approx(w_max, rel=tolerance)
    assert max(s.get_bending_moments()['values']) == pytest.approx(M_max, rel=tolerance)
    assert max(s.get_shear_forces()['values']) == pytest.approx(S_max, rel=tolerance)

def test_uniformly_distributed_load():
   
    q = 10
    
    M_max = q * L ** 2 / 8  # maximum moment
    S_max = q * L/2 # max shearing force
    w_max = -5 * q * L ** 4 / (384 * E * I)  # max displacement
    tolerance = 1e-4 #set  a tolerance of 0.01%

    s = Structure('test')

    s.add_beam(coord=[0, L],  E=E, I=I)
    s.add_distributed_load((q, q), (0, L))
    s.add_nodal_support({'uz': 0, 'ur': "NaN"}, 0)
    s.add_nodal_support({'uz': 0, 'ur': "NaN"}, L)
    s.add_nodes(200)
    s.add_elements(s.nodes)
    s.solve(s.build_global_matrix(), s.build_load_vector(), s.get_boudary_conditions())
  
    assert min(s.get_displacements()['vertical_displacements']) == pytest.approx(w_max, rel=tolerance)
    assert max(s.get_bending_moments()['values']) == pytest.approx(M_max, rel=tolerance)
    assert max(s.get_shear_forces()['values']) == pytest.approx(S_max, rel=1e-2)