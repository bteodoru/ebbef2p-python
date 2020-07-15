import pytest

from ebbef2p import Beam, Structure, NodalLoad, DistributedLoad, VlasovFoundationParameters


"""
Girija Vallabhan, C. V., & Das, Y. C. (1991). 
Modified Vlasov model for beams on elastic foundations. 
Journal of geotechnical engineering, 117(6), 956-966.
"""

L = 100
E = 432000
h = 3
w = 1
I = w * h**3 / 12
nu = 0.2
P = 20
q = 3



def test_concentrated_load_at_center():

    gamma = 0.507 #gamma parameter from reference paper
    tolerance = 1e-3 #set  a tolerance of 0.1%

    s = Structure('test')
        
    s.add_beam(Beam(coord=[0, L],  E=E, h=h, w=w))
    s.add_load(NodalLoad(P, L/2, 'fz'))
   
    s.discretize(50)

    vlasov_parameters = VlasovFoundationParameters({'E1': E/3000, 'E2': E/3000}, nu, L*0.25, gamma)
    p = s.compute_elastic_foundation_parameters({'E1': E/3000, 'E2': E/3000}, nu, L*0.25)

    assert p[-1].k == pytest.approx(vlasov_parameters.k, rel=tolerance)
    assert p[-1].t == pytest.approx(vlasov_parameters.t, rel=tolerance)

def test_concentrated_loads_at_ends():

    gamma = 1.239 #gamma parameter from reference paper
    tolerance = 1e-3 #set  a tolerance of 0.1%

    s = Structure('test')
        
    s.add_beam(Beam(coord=[0, L],  E=E, h=h, w=w))
    s.add_load(NodalLoad(P, 0, 'fz'))
    s.add_load(NodalLoad(P, L, 'fz'))
   
    s.discretize(50)

    vlasov_parameters = VlasovFoundationParameters({'E1': E/3000, 'E2': E/3000}, nu, L*0.25, gamma)
    p = s.compute_elastic_foundation_parameters({'E1': E/3000, 'E2': E/3000}, nu, L*0.25)

    assert p[-1].k == pytest.approx(vlasov_parameters.k, rel=tolerance)
    assert p[-1].t == pytest.approx(vlasov_parameters.t, rel=tolerance)    

def test_uniformly_distributed_load():

    gamma = 0.395 #gamma parameter from reference paper
    tolerance = 1e-3 #set  a tolerance of 0.1%

    s = Structure('test')
        
    s.add_beam(Beam(coord=[0, L],  E=E, h=h, w=w))
    s.add_load(DistributedLoad(value=(q, q), position=(0, L)))

    s.discretize(50)

    vlasov_parameters = VlasovFoundationParameters({'E1': E/3000, 'E2': E/3000}, nu, L*0.25, gamma)
    p = s.compute_elastic_foundation_parameters({'E1': E/3000, 'E2': E/3000}, nu, L*0.25)

    assert p[-1].k == pytest.approx(vlasov_parameters.k, rel=tolerance)
    assert p[-1].t == pytest.approx(vlasov_parameters.t, rel=0.015) #set  a tolerance of 1.5%