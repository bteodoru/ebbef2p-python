from ebbef2p.structure import Structure
from ebbef2p.vlasov_foundation_parameters import VlasovFoundationParameters

import numpy as np
import matplotlib.pyplot as plt
from ebbef2p.helpers import *

s = Structure('test')

s.add_beam(coord=[0, 10],  E=1, I=1)
s.add_nodal_load(50, 5, 'fz')
s.add_nodes(25)

#print(s.nodes)
p = s.compute_elastic_foundation_parameters({'E1': 10000, 'E2': 10000}, 0.3, 5)
#print(p)

for pa in p:
    pass
  # print(pa)


sol = s.solve(s.build_global_matrix(), s.build_load_vector(), s.get_boudary_conditions())

print(s.get_displacements()['vertical_displacements'][np.where(s.nodes == 5)[0][0]])
#print(np.where(s.nodes == 5)[0][0])