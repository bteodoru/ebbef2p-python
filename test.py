import pytest
import numpy as np
from ebbef2p.structure import Structure
from ebbef2p.nodal_support import NodalSupport


s = Structure('test')

s.add_beam(coord=[0, 25],  E=423000, I=2.25)
# s.add_nodal_load(P, 1, 'fz')
s.add_nodal_support({'uz': 0, 'ur': "NaN"})
#s.add_nodal_support(NodalSupport({'uz': 0, 'ur': "NaN"}, 0))
s.add_distributed_load((50, 50), (0, 5))
s.add_distributed_load((50, 60), (7, 12))
s.add_elastic_foundation([10, 0], [0, 25], 'k')
s.add_elastic_foundation([0, 0], [0, 25], 't')
s.add_nodes(50)
s.add_elements(s.nodes)
print(s.nodal_supports)
s.solve(s.build_global_matrix(), s.build_load_vector(), s.get_boudary_conditions())
