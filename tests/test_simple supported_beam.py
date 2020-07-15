import pytest
import numpy as np
from ebbef2p import Beam, Structure, NodalLoad, DistributedLoad, NodalSupport

L = 2
E = 1
h = 1
w = 1
I = w * h**3 / 12

def test_center_load():

    P = 100
       
    M_max = -P * L / 4  # maximum moment
    S_max = P/2 # max shearing force
    w_max = -P * L ** 3 / (48 * E * I)  # max displacement
    tolerance = 1e-6 #set  a tolerance of 0.0001%

    s = Structure('test')

    s.add_beam(Beam(coord=[0, L],  E=E, h=h, w=w))
    s.add_load(NodalLoad(P, L/2, 'fz'))
    s.add_nodal_support(NodalSupport(position=0, uz=0))
    s.add_nodal_support(NodalSupport(position=L, uz=0))
    s.discretize(25)
    s.add_elements()
    s.solve()
  
    assert min(s.get_displacements()['vertical_displacements']) == pytest.approx(w_max, rel=tolerance)
    assert min(s.get_bending_moments()['values']) == pytest.approx(M_max, rel=tolerance)
    assert max(s.get_shear_forces()['values']) == pytest.approx(S_max, rel=tolerance)

def test_uniformly_distributed_load():
   
    q = 10
    
    M_max = -q * L**2 / 8  # maximum moment
    S_max = q * L/2 # max shearing force
    w_max = -5 * q * L**4 / (384 * E * I)  # max displacement
    tolerance = 1e-4 #set  a tolerance of 0.01%

    s = Structure('test')

    s.add_beam(Beam(coord=[0, L],  E=E, h=h, w=w))
    s.add_load(DistributedLoad(value=(q, q), position=(0, L)))
    s.add_nodal_support(NodalSupport(position=0, uz=0))
    s.add_nodal_support(NodalSupport(position=L, uz=0))
    s.discretize(200)
    s.add_elements()
    s.solve()
  
    assert min(s.get_displacements()['vertical_displacements']) == pytest.approx(w_max, rel=tolerance)
    assert min(s.get_bending_moments()['values']) == pytest.approx(M_max, rel=tolerance)
    assert max(s.get_shear_forces()['values']) == pytest.approx(S_max, rel=1e-2)