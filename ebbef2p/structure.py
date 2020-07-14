from .beam import Beam
from .element import BeamElement
from .elastic_foundation import ElasticFoundation
from .vlasov_foundation_parameters import VlasovFoundationParameters
from .distributed_load import DistributedLoad
from .nodal_load import NodalLoad
from .nodal_support import NodalSupport
from .helpers import *
import numpy as np
from functools import reduce
import math
import matplotlib.pyplot as plt


class Structure():
    """Any model starts with an object of this class.

    Args:
        name (:obj:`str`): Name of the model
        beams (:obj:`list`): List of the beams 
        nodal_loads (:obj:`list`): List of the nodal loads
        distributed_loads (:obj:`list`): List of the distributed loads
        nodal_supports (:obj:`list`): List of the nodal constraints
        elastic_foundation (:obj:`list`): List 
        nodes (:obj:`list`): List of the nodes 
        elements (:obj:`list`): List of the finite elements 
    """

    def __init__(self, name):
        """Inits the Structure class."""

        self.name = name
        self.beams = []
        #self.forces = []
        #self.moments = []
        self.nodal_loads = []
        self.distributed_loads = []
        self.nodal_supports = []
        self.constraints = []
        #self.soil_conditions = []
        self.elastic_foundation = []
        self.nodes = []
        self.elements = []
        self.u = []
        self.forces = []
        #self.qn = []


    def get_key_points(self):
        """Key points are entities that geometrically marks the
        beams domains, nodal constraints, loads or elastic foundation.

        Returns:
            :obj:`list`: A list containing positions of the key points 
            for the current structure.
        """

        key_points = list(set().union(
            [val for sublist in 
            [b.coord for b in self.beams] for 
            val in sublist], 
            [nls.position for nls in self.nodal_loads 
            + self.nodal_supports],
            [val for sublist in 
            [q.position for q in self.distributed_loads] for 
            val in sublist],
            [val for sublist in 
            [ef.position for ef in self.elastic_foundation] for
            val in sublist])) 
        key_points.sort()
        return key_points

    def add_beam(self, beam):
        """Add a beam to the structure model.

        Args:
            beam (:class:`~ebbef2p.beam.Beam`): Beam instance.

        Raises:
            AttributeError: If ``beam`` is not an instance of 
                :class:`~ebbef2p.beam.Beam`.
        """

        if isinstance(beam, Beam):
            self.beams.append(beam)
        else:
            raise AttributeError('The given parameter is not an ' 
                                 'instance of Beam')


        #self.beams.append(Beam(coord, E, I))

    def add_load(self, load):
        """Add a load to the structure model.

        Args:
            load (:class:`~ebbef2p.nodal_load.NodalLoad` or \
            :class:`~ebbef2p.distributed_load.DistributedLoad`):
                Nodal load or Distributed load instance.
        
        Raises:
            AttributeError: If ``load`` is not an instance of 
                either :class:`~ebbef2p.nodal_load.NodalLoad`
                or :class:`~ebbef2p.distributed_load.DistributedLoad`
        """

        if isinstance(load, NodalLoad):
            self.nodal_loads.append(load)
        elif isinstance(load, DistributedLoad):
            self.distributed_loads.append(load)
        else:
            raise AttributeError('The given parameter is not an ' 
                                 'instance of NodalLoad or '
                                 'DistributedLoad')

    def add_nodal_support(self, support):
        """Add nodal support to the structure model.

        Args:
            support (:class:`~ebbef2p.nodal_support.NodalSupport`): 
                Nodal support instance.

        Raises:
            AttributeError: If ``support`` is not an instance of 
                :class:`~ebbef2p.nodal_support.NodalSupport`
        """   

        if isinstance(support, NodalSupport):
            self.nodal_supports.append(support)
        else:
            raise AttributeError('The given parameter is not an '
                                  'instance of NodalSupport')
  
    def add_elastic_foundation(self, elastic_foundation):
        """Add elastic foundation support to the structure model.

        Args:
            elastic_foundation (:class:`~ebbef2p.elastic_foundation.ElasticFoundation`):
                Elastic foundation instance.

        Raises:
            AttributeError: If ``elastic_foundation`` is not an
                instance of 
                :class:`~ebbef2p.elastic_foundation.ElasticFoundation`
        """

        if isinstance(elastic_foundation, ElasticFoundation):   
            self.elastic_foundation.append(elastic_foundation)
        else:
            raise AttributeError('The given parameter is not an '
                                  'instance of ElasticFoundation')
        
    def discretize(self, n):
        """Takes an already defined 
        :class:`~ebbef2p.structure.Structure` object and divide into at 
        least ``n`` elements.

        Args:
            n (:obj:`int`): Minimum number of elements to divide the 
                structure.
        """

        delta = (max(self.get_key_points()) - min(self.get_key_points())) / n
        nodes = np.array([])
        for i in range(1, len(self.get_key_points())):
            nodes = np.union1d(nodes, np.linspace(
                self.get_key_points()[i-1], self.get_key_points()[i], 
                1 + math.ceil((self.get_key_points()[i]
                - self.get_key_points()[i-1]) / delta)))
        self.nodes = nodes

    def add_elements(self):
        """Add beam finite elements to the current structure model
        """
        elements = []
        for x in pairwise(self.nodes):
            ni, nj = x[0], x[1]
            ki, kj = 0, 0
            ti, tj = 0, 0
            for b in self.beams:
                if is_within(x, b.coord):
                    E = b.E
                    h = b.h
                    w = b.w
            for ef in self.elastic_foundation:
                if ef.type == 'k':
                    if is_within(x, ef.position):
                        ki += np.interp(ni, ef.position, ef.value)
                        kj += np.interp(nj, ef.position, ef.value)
                if ef.type == 't':
                    if is_within(x, ef.position):
                        ti += np.interp(ni, ef.position, ef.value)
                        tj += np.interp(nj, ef.position, ef.value)

            beam_element = BeamElement((ni, nj), E, h, w, [ki * w, kj * w], [ti * w, tj * w])
            elements.append(beam_element)

        self.elements = elements

    def compute_elastic_foundation_parameters(self, Edef, nu, depth):
        """Compute elastic foundation parameters

        Args:
            Edef (:obj:`float`): Modulus of elasticity
            nu (:obj:`float`): Poisson's ratio
            depth (:obj:`float`): Foundation depth

        Returns:
            :obj:`list`: List of 
            :class:`~ebbef2p.vlasov_foundation_parameters.VlasovFoundationParameters`
            instances.
        """

        gamma = 1
        gamma_it = []
        parameters = []
        structure_length = sum([b.length for b in self.beams])
        while gamma > 0:
            vlasov_parameters = VlasovFoundationParameters(
                    Edef, nu, depth, gamma)
            gamma_it.append(gamma)
            parameters.append(vlasov_parameters)
            self.elastic_foundation = []
            self.add_elastic_foundation(ElasticFoundation(
                    [vlasov_parameters.k, vlasov_parameters.k], 
                    [0, structure_length], 'k'))
            self.add_elastic_foundation(ElasticFoundation(
                    [vlasov_parameters.t, vlasov_parameters.t], 
                    [0, structure_length], 't'))
             
            self.add_elements()
            self.solve()


            gamma = vlasov_parameters.get_gamma(
                self.get_vertical_displacements(), self.nodes, 
                vlasov_parameters.k, vlasov_parameters.t)
            if len(gamma_it) > 2:
                    if abs(gamma_it[-2]-gamma_it[-1]) < 0.0001:
                        break
        
        return parameters

    def build_global_matrix(self):
        """Assembles the global stiffness

        Returns:
            :obj:`numpy.array`: The global stiffness matrix
        """
        elements = self.elements
        K = np.zeros((2*(len(elements)+1), 2*(len(elements)+1)))
        for e, element in enumerate(elements):
            K[e*2:e*2+4, e*2:e*2+4] += element.stiffness

        return K

    def get_forces(self):
        f = np.zeros((len(self.elements), 4))
        for e, element in enumerate(self.elements):
            #f = np.append(f, np.dot(element.ke, self.u[2*e:2*e+4]))
            # print(f)
            f[e] = element.forces(self.u[2*e:2*e+4].flatten())

        return f

        #self.forces = f

    def get_vertical_displacements(self):
        """Return the nodal displacements based on the analysis results.

        Returns:
            :obj:`dict`: A dictionary of {:obj:`str`: :obj:`list`} where the 
            keys are ```vertical_displacements``` and  ```rotations```.       
        """

        return dict(vertical_displacements=[-w for w in self.u.flatten()[0::2]], rotations=[-t for t in self.u.flatten()[1::2]])

    def get_shear_forces(self):
        coords = list(itertools.chain(*[e.coord for e in self.elements]))
        shear_forces = list(itertools.chain(
           # *[(f[0], -f[2]) for f in self.get_forces()]))
            *[(f[0], -f[2]) for f in self.get_forces()]))
        # return [coords, shear_forces]
        return dict(coords=coords, values=shear_forces)
        
    def get_bending_moments(self):
        coords = list(itertools.chain(*[e.coord for e in self.elements]))
        bending_moments = list(itertools.chain(
            #*[(-f[1], f[3]) for f in self.get_forces()]))
            *[(-f[1], f[3]) for f in self.get_forces()]))

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

    def get_boudary_conditions(self):
        bc = np.array([])
        # constraints: {ur: 0, uz: 0}, position: 0 ##{'uz': 0, 'ur': "null"}, 9
        for c in self.nodal_supports:
            ni = np.where(self.nodes == c.position)[0]
            #print(bool(c.constraints.get('uz') != 'NaN'))

            if c.uz is not None:
                # if 'uz' in c.constraints:
                bc = np.append(bc, [ni[0]*2, c.uz])
            if c.ur is not None:
                # if 'ur' in c.constraints:
                bc = np.append(bc, [1+ni[0]*2, c.ur])
        return np.reshape(bc, (int(0.5*len(bc)), 2))
        # print(bc)

    def build_load_vector(self):
        """Build load vector 

        Returns:
            :obj:`numpy array`: The load vector
        """
        lv = np.zeros(len(self.nodes)*2)

        for nl in self.nodal_loads:
            index = np.where(self.nodes == nl.position)[0]
            # print(index)
            if nl.type == 'fz':
                lv[2*index] = nl.value
            if nl.type == 'my':
                lv[2*index + 1] = nl.value

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

    def solve(self):
        """Solve static FE-equations considering boundary conditions.
        Adapted from CALFEM (for Python): 
        https://github.com/CALFEM/calfem-python

        Returns:
            :obj:`tuple`: solution including boundary values (``a``) and 
            reaction force vector (``Q``), 
            dim(a)=dim(Q)= nd x 1, nd : number of dof's
        """
        # """
        # Solve static FE-equations considering boundary conditions.
        # Adapted from CALFEM (for Python): 
        # https://github.com/CALFEM/calfem-python

        # Parameters:

        # K            global stiffness matrix, dim(K)= nd x nd
        # f            global load vector, dim(f)= nd x 1

        # constraints TO DO
        # Returns:

        # a           solution including boundary values
        # Q           reaction force vector
        #             dim(a)=dim(Q)= nd x 1, nd : number of dof's

        # """

        K = self.build_global_matrix()
        f = self.build_load_vector()
        constraints = self.get_boudary_conditions()
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


        Q = K*np.asmatrix(a)-f.reshape(nDofs, 1)

        self.u = a

        return (np.asmatrix(a), Q)

       

    def __str__(self):
        """[summary]

        Returns:
            [type]: [description]
        """
        return f"Beams: {self.beams} \
                \nNodal Loads: {self.nodal_loads} \
                \nDistributed Loads: {self.distributed_loads} \
                \nElastic Foundation: {self.elastic_foundation}"
