import pytest

from ebbef2p.structure import Structure
from ebbef2p.vlasov_foundation_parameters import VlasovFoundationParameters

def test_elastic_foundation():
    s = Structure('test')

    s.add_beam(coord=[0, 10],  E=1, I=1)
    s.add_nodal_load(50, 5, 'fz')
   
    s.add_nodes(1)
    s.add_elements(s.nodes)
    s.compute_elastic_foundation_parameters({'E1': 10000, 'E2': 10000}, 0.3, 5)

    K=s.build_global_matrix()
    sol = s.solve(K, s.build_load_vector(), s.get_boudary_conditions())
    #print(s.get_shear_forces()['values'][0])
    print(s.get_displacements()['vertical_displacements'][0])
   
    
    
    #assert s.get_shear_forces()['values'][0] == pytest.approx(25, 0.1)
    #assert s.get_bending_moments()['values'][1] == pytest.approx(125, 0.1)
    #assert s.get_displacements()['vertical_displacements'][0] == pytest.approx(0, 0.1)