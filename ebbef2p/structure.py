from .beam import Beam
from .element import BeamElement
from .soil_conditions import SoilCondition
from .distributed_load import DistributedLoad
from .nodal_load import NodalLoad
from .nodal_support import NodalSupport
from .helpers import *
import numpy as np
from functools import reduce
import math
import matplotlib.pyplot as plt


class Structure():
    """docstring for Structure"""

    def __init__(self, name):
        """Constructor method
        """        
        self.name = name
        self.beams = []
        #self.forces = []
        #self.moments = []
        self.nodal_loads = []
        self.distributed_loads = []
        self.nodal_supports = []
        self.constraints = []
        self.soil_conditions = []
        self.nodes = []
        self.u = []
        self.forces = []
        #self.qn = []

    # @property
    def get_points_(self):
                
        v = np.array([])
        for b in self.beams:
            # [2, 3, 5] -> [0, 2, 5, 10] np.cumsum()
            v = np.union1d(v, b.coord)
        # for f in self.forces:
        #     v = np.union1d(v, f[1])
        for nl in self.nodal_loads:
            v = np.union1d(v, nl.position)

        for ns in self.nodal_supports:
            v = np.union1d(v, ns.position)
        # for m in self.moments:
        #     v = np.union1d(v, m[1])
        for q in self.distributed_loads:
            v = np.union1d(v, q.position)
        for sc in self.soil_conditions:
            #v = reduce(np.union1d, (v, sc.k[1], sc.t[1]))
            #reduce(np.union1d, ([1, 3, 4, 3], [3, 1, 2, 1], [6, 3, 4, 2]))
            v = np.union1d(v, sc.position)

        return v

    def get_points(self):
        v = list(set().union([val for sublist in [b.coord for b in self.beams] for val in sublist],
                             [nls.position for nls in self.nodal_loads +
                                 self.nodal_supports],
                             [val for sublist in [
                                 q.position for q in self.distributed_loads] for val in sublist],
                             [val for sublist in [sc.position for sc in self.soil_conditions] for val in sublist])) \
            # [val for sublist in [sc.t[1] for sc in self.soil_conditions] for val in sublist] ))
        v.sort()
        return v

    def add_beam(self, coord, E, I):
        self.beams.append(Beam(coord, E, I))

    def add_nodal_load(self, value, position, type):
        self.nodal_loads.append(NodalLoad(value, position, type))

    def add_distributed_load(self, value, position):
        self.distributed_loads.append(DistributedLoad(value, position))

    def add_nodal_support(self, constraints, position):
        self.nodal_supports.append(NodalSupport(constraints, position))

    def add_soil_condition(self, value, position, type):
        # self.soil_conditions[0].append(k)
        # self.soil_conditions[1].append(t)
        self.soil_conditions.append(SoilCondition(value, position, type))

    def discretization(self, n):
        points = self.get_points()
        size = (max(points) - min(points))/n
        nodes = np.array([])
        #nodes = []

        for i in range(1, len(points)):
            nodes = np.union1d(nodes, np.linspace(
                points[i-1], points[i], 1 + math.ceil((points[i]-points[i-1])/size)))
            #nodes.append(np.linspace(points[i-1], points[i], 1 + math.ceil((points[i]-points[i-1])/size)))
        # return nodes
        self.nodes = nodes
        self.elements = self.add_elements(nodes)
        return nodes

    def mesh(self, n):
        points = self.get_points()
        size = (max(points) - min(points))/n
        nodes = np.array([])
        #nodes = []

        for i in range(1, len(points)):
            nodes = np.union1d(nodes, np.linspace(
                points[i-1], points[i], 1 + math.ceil((points[i]-points[i-1])/size)))
            #nodes.append(np.linspace(points[i-1], points[i], 1 + math.ceil((points[i]-points[i-1])/size)))
        # return nodes
        self.nodes = nodes
        #self.elements = self.add_elements(nodes)
        return nodes

    def add_elements(self, nodes):
        elements = []
        for x in pairwise(nodes):
            ni, nj = x[0], x[1]
            for b in self.beams:
               # if b.coord[0] <  nj <= b.coord[1]:
                if is_within(x, b.coord):
                    E = b.E
                   # print(E)
                    I = b.I
            for sc in self.soil_conditions:
                # print(sc)
                # if sc.k[1][0] < nj <= sc.k[1][1]:
                if sc.type == 'k':
                    if is_within(x, sc.position):
                        ki = np.interp(ni, sc.position, sc.value)
                        kj = np.interp(nj, sc.position, sc.value)
                    else:
                        ki = 0
                        kj = 0
                if sc.type == 't':
                    if is_within(x, sc.position):
                        ti = np.interp(ni, sc.position, sc.value)
                        tj = np.interp(nj, sc.position, sc.value)
                    else:
                        ti = 0
                        tj = 0

                # if is_within(x, sc.position):
                #     if sc['type'] == 'k':
                #         ki = np.interp(ni, sc['position'], sc['value'])
                #         kj = np.interp(nj, sc['position'], sc['value'])
                #     if sc['type'] == 't':
                #         ti = np.interp(ni, sc['position'], sc['value'])
                #         tj = np.interp(nj, sc['position'], sc['value'])
                # if is_within(x, sc.k[1]):

                #     ki = np.polyfit(sc.k[1], sc.k[0], 1)[
                #         0]*ni + np.polyfit(sc.k[1], sc.k[0], 1)[1]
                #     kj = np.polyfit(sc.k[1], sc.k[0], 1)[
                #         0]*nj + np.polyfit(sc.k[1], sc.k[0], 1)[1]
                # else:
                #     ki = 0
                #     kj = 0

                # ## if sc.t[1][0] < nj <= sc.t[1][1]:
                # if is_within(x, sc.t[1]):
                #     ti = np.polyfit(sc.t[1], sc.t[0], 1)[
                #         0]*ni + np.polyfit(sc.t[1], sc.t[0], 1)[1]
                #     tj = np.polyfit(sc.t[1], sc.t[0], 1)[
                #         0]*nj + np.polyfit(sc.t[1], sc.t[0], 1)[1]
                # else:
                #     ti = 0
                #     tj = 0

            be = BeamElement((ni, nj), E, I, [ki, kj], [ti, tj])
            elements.append(be)

        return elements

    def build_global_matrix(self):
        elements = self.elements
        K = np.zeros((2*(len(elements)+1), 2*(len(elements)+1)))
        for e, element in enumerate(elements):
            # ke = flexural_stiffness(element)
            K[e*2:e*2+4, e*2:e*2+4] += element.stiffness

        return K

    def build_global_flexural_matrix(self, elements):
        K = np.zeros((2*(len(elements)+1), 2*(len(elements)+1)))
        for e, element in enumerate(elements):
            # ke = flexural_stiffness(element)
            K[e*2:e*2+4, e*2:e*2+4] += element.ke

        return K

    def get_forces(self):
        f = np.zeros((len(self.elements), 4))
        for e, element in enumerate(self.elements):
            #f = np.append(f, np.dot(element.ke, self.u[2*e:2*e+4]))
            # print(f)
            f[e] = element.forces(self.u[2*e:2*e+4].flatten())

        return f

        #self.forces = f

    def get_displacements(self):

        #return dict(vertical_displacements=[w for w in self.u.flatten()[0::2]], rotations=[t for t in self.u.flatten()[1::2]])
        return dict(vertical_displacements=[-w for w in self.u.flatten()[0::2]], rotations=[-t for t in self.u.flatten()[1::2]])

    def get_shear_forces(self):
        coords = list(itertools.chain(*[e.coord for e in self.elements]))
        shear_forces = list(itertools.chain(
           # *[(f[0], -f[2]) for f in self.get_forces()]))
            *[(-f[0], f[2]) for f in self.get_forces()]))
        # return [coords, shear_forces]
        return dict(coords=coords, values=shear_forces)
        

    def get_bending_moments(self):
        coords = list(itertools.chain(*[e.coord for e in self.elements]))
        bending_moments = list(itertools.chain(
            #*[(-f[1], f[3]) for f in self.get_forces()]))
            *[(f[1], -f[3]) for f in self.get_forces()]))

        return dict(coords=coords, values=bending_moments)

    def plot_shear_diagram(self):
        #coords = list(itertools.chain(*[e.coord for e in self.elements]))
        #shearing_forces = list(itertools.chain(*[(f[0], -f[2]) for f in self.get_forces()]))
       # print(shearing_forces)

        plt.plot(self.get_shear_forces()['coords'],
                 self.get_shear_forces()['values'])
        plt.show()

    def plot_moment_diagram(self):
        #coords = list(itertools.chain(*[e.coord for e in self.elements]))
        #bending_moments = list(itertools.chain(*[(-f[1], f[3]) for f in self.get_forces()]))

        plt.plot(self.get_bending_moments()[
                 'coords'], self.get_bending_moments()['values'])
        plt.show()

    def get_bc(self):
        #cs = np.zeros(len(self.nodes)*2)
        bcdof = np.array([])
        for c in self.constraints:
           # ni = np.where(self.nodes == c[2])[0]
            ni = np.where(self.nodes == c['x'])[0]
            # print(ni)
            if 'uz' in c:
               # bcdof.append([int(ni[0]*2), c['ux']])
               #{'uz': 0, 'x': 9}
                bcdof = np.append(bcdof, [ni[0]*2, c['uz']])
            if 'ur' in c:
                bcdof = np.append(bcdof, [1+ni[0]*2, c['ur']])
               # bcdof.append([int(1+ni[0]*2), c['ur']])
           # bcdof.append([ni[0], c[0], c[1]])
           # bcdof.append([ni[0]+1, c[0], c[1]])
            # print(index)
            #bcdof[2*index] = c[0]
            #bcdof[2*index + 1] = c[1]

        return np.reshape(bcdof, (len(self.constraints), 2))

    def get_boudary_conditions(self):
        bc = np.array([])
        # constraints: {ur: 0, uz: 0}, position: 0 ##{'uz': 0, 'ur': "null"}, 9
        for c in self.nodal_supports:
            ni = np.where(self.nodes == c.position)[0]
            #print(bool(c.constraints.get('uz') != 'NaN'))

            if c.constraints['uz'] != 'NaN':
                # if 'uz' in c.constraints:
                bc = np.append(bc, [ni[0]*2, c.constraints['uz']])
            if c.constraints['ur'] != 'NaN':
                # if 'ur' in c.constraints:
                bc = np.append(bc, [1+ni[0]*2, c.constraints['ur']])
        return np.reshape(bc, (int(0.5*len(bc)), 2))
        # print(bc)

    def build_displacements_vector(self):
        pass

    def build_load_vector(self):
        lv = np.zeros(len(self.nodes)*2)
        # for f in self.forces:
        #     index = np.where(self.nodes == f[1])[0]
        #     # print(2*index)
        #     lv[2*index] = f[0]
        # for m in self.moments:
        #     index = np.where(self.nodes == m[1])[0]
        #     lv[2*index + 1] = m[0]

        # for i, n in enumerate(self.nodes):
        #     for q in self.distributed_loads:
        #         if q.position[0] <= n <= q.position[1]:
        # qn = np.zeros(len(self.nodes))
        #             qn[i] = np.interp(n, q.position, q.value)
        #lv[2*i] += (n[i+1] - n[i])/20*(7*qn[i] + 3*qn[i+1])
        #self.qn = qn

        for nl in self.nodal_loads:
            index = np.where(self.nodes == nl.position)[0]
            # print(index)
            if nl.type == 'fz':
                lv[2*index] = nl.value
            if nl.type == 'my':
                lv[2*index + 1] = nl.value

        #t = np.zeros(len(self.nodes)*2)
        # for ni, q in enumerate(pairwise(qn)):

        #     qi, qj = q[0], q[1]
        #     l = self.nodes[ni+1] - self.nodes[ni]
        #     F1q = l/20 * (7*qi + 3*qj)
        #     M1q = l**2/60 * (3*qi + 2*qj)
        #     F2q = l/20 * (3*qi + 7*qj)
        #     M2q = -l**2/60 * (2*qi + 3*qj)

            #t += [F1q, M1q, F2q, M2q]
        for e, n in enumerate(pairwise(self.nodes)):
            for q in self.distributed_loads:
                if is_within(n, q.position):
                    qi = np.interp(n[0], q.position, q.value)
                    qj = np.interp(n[1], q.position, q.value)
                    #print(f"qi: {qi} \nqj: {qj}")
                    #l = (nj - ni)
                    F1q = (n[1] - n[0])/20 * (7*qi + 3*qj)
                    M1q = (n[1] - n[0])**2/60 * (3*qi + 2*qj)
                    F2q = (n[1] - n[0])/20 * (3*qi + 7*qj)
                    M2q = -(n[1] - n[0])**2/60 * (2*qi + 3*qj)

                    lv[2*e:2*e+4] += [F1q, M1q, F2q, M2q]

        return lv

    def solve(self, K, f, constraints):
        """
        Solve static FE-equations considering boundary conditions.
        Adapted from CALFEM (for Python): https://github.com/CALFEM/calfem-python

        Parameters:

        K            global stiffness matrix, dim(K)= nd x nd
        f            global load vector, dim(f)= nd x 1

        constraints TO DO
        Returns:

        a           solution including boundary values
        Q           reaction force vector
                    dim(a)=dim(Q)= nd x 1, nd : number of dof's

        """

        nDofs = K.shape[0]
        bcVal = constraints[:, 1]
        bcPrescr = constraints[:, 0].astype(int)
        nPdofs = bcPrescr.shape[0]

        bc = np.ones(nDofs, 'bool')
        bcDofs = np.arange(nDofs)

        bc[np.ix_(bcPrescr)] = False
        bcDofs = bcDofs[bc]

        fsys = f[bcDofs]-K[np.ix_((bcDofs), (bcPrescr))].dot(bcVal)
        asys = np.linalg.solve(K[np.ix_((bcDofs), (bcDofs))], fsys)

        a = np.zeros([nDofs, 1])
        a[np.ix_(bcPrescr)] = np.asmatrix(bcVal).reshape(nPdofs, 1)
        a[np.ix_(bcDofs)] = np.asmatrix(asys).reshape(bcDofs.shape[0], 1)
        #print(f"K.shape = {asys}")
        #print(f"a shape = {a.shape}")
        #print(f"f shape = {f.shape}")
        # print(a.flatten())

        Q = K*np.asmatrix(a)-f.reshape(nDofs, 1)

        self.u = a

        return (np.asmatrix(a), Q)

       
    def elements_forces(self):
        f = np.zeros(len(self.nodes)*2)
        for i, el in enumerate(self.elements):
            f[2*i:2*i+4] += el.forces(self.u[2*i:2*i+4].flatten())

        self.forces = f

    def __str__(self):
        return f"Beams: {self.beams} \nNodal Loads: {self.nodal_loads} \nDistributed Loads: {self.distributed_loads} \nSoil Conditions: {self.soil_conditions}"
