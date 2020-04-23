import pytest

from ebbef2p.structure import Structure

def test_simple_supported_beam_P():
    s = Structure('test')

    s.add_beam(coord=[0, 10],  E=1, I=1)
    s.add_nodal_load(50, 5, 'fz')
    s.add_nodal_support({'uz': 0, 'ur': "NaN"}, 0)
    s.add_nodal_support({'uz': 0, 'ur': "NaN"}, 10)
    s.add_soil_condition([0, 0], [0, 10], 'k')
    s.add_soil_condition([0, 0], [0, 10], 't')
    s.discretization(1)
    K=s.build_global_matrix()
    sol = s.solve(K, s.build_load_vector(), s.get_boudary_conditions())
    #print(s.get_shear_forces()['values'][0])
   
    
    
    assert s.get_shear_forces()['values'][0] == pytest.approx(25, 0.1)
    assert s.get_bending_moments()['values'][1] == pytest.approx(125, 0.1)

def test_simple_supported_beam_q():
    s = Structure('test')

    s.add_beam(coord=[0, 10],  E=1, I=1)
    s.add_distributed_load((50, 50), (0, 10))
    s.add_nodal_support({'uz': 0, 'ur': "NaN"}, 0)
    s.add_nodal_support({'uz': 0, 'ur': "NaN"}, 10)
    s.add_soil_condition([0, 0], [0, 10], 'k')
    s.add_soil_condition([0, 0], [0, 10], 't')
    s.discretization(250)
    K=s.build_global_matrix()
    sol = s.solve(K, s.build_load_vector(), s.get_boudary_conditions())
  
   
    
    
    assert s.get_shear_forces()['values'][0] == pytest.approx(250, 0.1)
    #assert s.get_bending_moments()['values'][1] == pytest.approx(125, 0.1)