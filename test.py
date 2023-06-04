import pytest
import numpy as np
from ebbef2p import Beam
from ebbef2p import NodalLoad
from ebbef2p import Structure
from ebbef2p import ElasticFoundation
from ebbef2p.nodal_support import NodalSupport

L = 100
E = 432000
I = 2.25
nu = 0.2
P = 20
q = 3


s = Structure('test')
#print(s.beams[0])
    
s.add_beam(Beam(coord=[0, 20],  E=27000000, w=0.5, h=1))
print(s.beams[0])
s.add_load(NodalLoad(1000, 10, 'fz'))
s.add_elastic_foundation(ElasticFoundation([4815, 4815], [0, 20], 'k'))     
s.add_elastic_foundation(ElasticFoundation([12676, 12676], [0, 20], 't'))               
print(f"c: {s.get_boudary_conditions()}")

s.discretize(25)
s.add_elements()
s.solve()
s.plot_moment_diagram()
p = s.compute_elastic_foundation_parameters({'E1': 20000, 'E2': 20000}, 0.25, 5)
print(p[-1])
#print(s.get_bending_moments()['values'])
#s.plot_moment_diagram()


s = Structure('test')

#s.add_beam(coord=[0, 25],  E=423000, I=2.25)
b = Beam(coord=[0, 20],  E=27000000, w=0.5, h=1)
nl = NodalLoad(value=1000, position=10, type='fz')
ns = NodalSupport(0, 0)
s.add_beam(b)
s.add_load(nl)
#s.add_nodal_support(ns)
#print(s)
#print(s.beams[0])
s.discretize(n=10)
s.add_elements()
sol = s.solve()
print(sol)
print(s.get_forces())
#print(s.get_key_points())
# s.add_nodal_load(P, 1, 'fz')
#s.add_nodal_support({'uz': 0, 'ur': "NaN"})
#s.add_nodal_support(NodalSupport({'uz': 0, 'ur': "NaN"}, 0))
#s.add_distributed_load((50, 50), (0, 5))
#s.add_distributed_load((50, 60), (7, 12))
#s.add_elastic_foundation([10, 0], [0, 25], 'k')
#s.add_elastic_foundation([0, 0], [0, 25], 't')
#s.add_nodes(50)
#s.add_elements(s.nodes)
#print(s.nodal_supports)
#s.solve(s.build_global_matrix(), s.build_load_vector(), s.get_boudary_conditions())
