from ebbef2p.structure import Structure
from ebbef2p.vlasov_foundation_parameters import VlasovFoundationParameters

import numpy as np
import matplotlib.pyplot as plt
from ebbef2p.helpers import *

s = Structure('test')

s.add_beam(coord=[0, 100],  E=432000, I=2.25)
#s.add_nodal_load(30, 50, 'fz')
#s.add_nodal_load(20, 100, 'fz')
s.add_distributed_load((3, 3), (0, 100))
s.add_nodes(50)

#print(s.nodes)
p = s.compute_elastic_foundation_parameters({'E1': 144, 'E2': 144}, 0.2, 25)
#print(p)

for pa in p:
    #pass
   print(pa)

#print(p[-1].k)
sol = s.solve(s.build_global_matrix(), s.build_load_vector(), s.get_boudary_conditions())
print(sol)
print(s.u)

print(s.get_displacements()['vertical_displacements'][np.where(s.nodes == 50)[0][0]])
print(s.get_displacements()['vertical_displacements'][0])
#print(np.where(s.nodes == 5)[0][0])