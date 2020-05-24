import pytest

from ebbef2p.structure import Structure
from ebbef2p.vlasov_foundation_parameters import VlasovFoundationParameters

"""
Girija Vallabhan, C. V., & Das, Y. C. (1991). 
Modified Vlasov model for beams on elastic foundations. 
Journal of geotechnical engineering, 117(6), 956-966.
"""

L = 100
E = 432000
I = 2.25
nu = 0.2
P = 20
q = 3



def test_concentrated_load_at_center():

    gamma = 0.507 #gamma parameter from reference paper
    tolerance = 1e-3 #set  a tolerance of 0.1%

    s = Structure('test')
        
    s.add_beam(coord=[0, L],  E=E, I=I)
    s.add_nodal_load(P, L/2, 'fz')
   
    s.add_nodes(50)

    vlasov_parameters = VlasovFoundationParameters({'E1': E/3000, 'E2': E/3000}, nu, L*0.25, gamma)
    p = s.compute_elastic_foundation_parameters({'E1': E/3000, 'E2': E/3000}, nu, L*0.25)

    assert p[-1].k == pytest.approx(vlasov_parameters.k, rel=tolerance)
    assert p[-1].t == pytest.approx(vlasov_parameters.t, rel=tolerance)

def test_concentrated_loads_at_ends():

    gamma = 1.239 #gamma parameter from reference paper
    tolerance = 1e-3 #set  a tolerance of 0.1%

    s = Structure('test')
        
    s.add_beam(coord=[0, L],  E=E, I=I)
    s.add_nodal_load(P, 0, 'fz')
    s.add_nodal_load(P, L, 'fz')
   
    s.add_nodes(50)

    vlasov_parameters = VlasovFoundationParameters({'E1': E/3000, 'E2': E/3000}, nu, L*0.25, gamma)
    p = s.compute_elastic_foundation_parameters({'E1': E/3000, 'E2': E/3000}, nu, L*0.25)

    assert p[-1].k == pytest.approx(vlasov_parameters.k, rel=tolerance)
    assert p[-1].t == pytest.approx(vlasov_parameters.t, rel=tolerance)    

def test_uniformly_distributed_load():

    gamma = 0.395 #gamma parameter from reference paper
    tolerance = 1e-3 #set  a tolerance of 0.1%

    s = Structure('test')
        
    s.add_beam(coord=[0, L],  E=E, I=I)
    s.add_distributed_load((q, q), (0, L))

    s.add_nodes(50)

    vlasov_parameters = VlasovFoundationParameters({'E1': E/3000, 'E2': E/3000}, nu, L*0.25, gamma)
    p = s.compute_elastic_foundation_parameters({'E1': E/3000, 'E2': E/3000}, nu, L*0.25)

    assert p[-1].k == pytest.approx(vlasov_parameters.k, rel=tolerance)
    assert p[-1].t == pytest.approx(vlasov_parameters.t, rel=0.015) #set  a tolerance of 1.5%