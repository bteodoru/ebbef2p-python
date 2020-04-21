from .beam import Beam
import numpy as np

class BeamElement(Beam):
    def __init__(self, coord, E, I, k, t):
        Beam.__init__(self, coord, E, I)
       # self.coord = coord
        #self.E = E
        #self.I = I
        self.k = k
        self.t = t
        self.ke = self.flexural_stiffness()
        self.k1 = self.first_parameter_foundation_stiffness()
        self.k2 = self.second_parameter_foundation_stiffness()
        self.stiffness = self.stiffness()
        #self.forces = self.forces()
   # def __init__(self):
        #self.coord = coord
   # @property
   # def length(self):
    #    return self

    def flexural_stiffness(self):
        E = self.E
        I = self.I
        l = self.length
        ke = (I*E/l**3)*np.array([
            [12,    6*l,  -12,    6*l],
            [6*l, 4*l**2, -6*l, 2*l**2],
            [-12,   -6*l,   12,   -6*l],
            [6*l, 2*l**2, -6*l, 4*l**2]])

        return ke

    def first_parameter_foundation_stiffness(self):
        k = self.k
        l = self.length
        k1 = np.array([
            [1/35*l*(10*k[0]+3*k[1]),      1/420*l**2*(15*k[0]+7*k[1]),
             9/140*l*(k[0]+k[1]),        -1/420*l**2*(7*k[0]+6*k[1])],
            [1/420*l**2*(15*k[0]+7*k[1]),  1/840*l**3*(5*k[0]+3*k[1]),
             1/420*l**2*(6*k[0]+7*k[1]), -1/280*l**3*(k[0]+k[1])],
            [9/140*l*(k[0]+k[1]),          1/420*l**2*(6*k[0]+7*k[1]),
             1/35*l*(3*k[0]+10*k[1]),    -1/420*l**2*(7*k[0]+15*k[1])],
            [-1/420*l**2*(7*k[0]+6*k[1]), -1/280*l**3*(k[0]+k[1]),     -1/420*l**2*(7*k[0]+15*k[1]), 1/840*l**3*(3*k[0]+5*k[1])]])

        return k1

    def second_parameter_foundation_stiffness(self):
        t = self.t
        l = self.length
        k2 = np.array([
            [3/5/l*(t[0]+t[1]),  1/10*t[1],            -
             3/5/l*(t[0]+t[1]),  1/10*t[0]],
            [1/10*t[1],           1/30*l *
             (3*t[0]+t[1]), -1/10*t[1],         -1/60*l*(t[0]+t[1])],
            [-3/5/l*(t[0]+t[1]), -1/10*t[1],
             3/5/l*(t[0]+t[1]), -1/10*t[0]],
            [1/10*t[0],          -1/60*l*(t[0]+t[1]),   -1/10*t[0],          1/30*l*(t[0]+3*t[1])]])

        return k2

    def stiffness(self):
        return self.ke + self.k1 + self.k2

    def is_within(self, rcoord):
        if self.coord[0] >= rcoord[0] and self.coord[1] <= rcoord[1]:
            return True
        else:
            return False

    def forces(self, u):
        return np.dot(self.ke, u)

    def __str__(self):
        return f"Coords: {self.coord} \nYoung's modulus: {self.E} \nLength: {self.length} \nk: {self.k} \nt: {self.t} \nforces: {self.forces}"
