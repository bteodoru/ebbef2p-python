import numpy as np

class VlasovFoundationParameters():

    def __init__(self, Edef, nu, depth, gamma):
        self.Edef = Edef
        self.nu = nu
        self.depth = depth
        self.gamma = gamma
        #self.k = None
        #self.t = None

    @property
    def k(self):
        return (((1 - self.nu) / ((1 + self.nu) * (1 - 2*self.nu) * self.depth))
                * (self.Edef['E1'] * (2*self.gamma * np.sinh(2*self.gamma)
                + 4 * (self.gamma**2)) + (self.Edef['E2'] - self.Edef['E1'])
                * (np.cosh(2*self.gamma) - 1 + 2*self.gamma**2)) 
                / (8 * ((np.sinh(self.gamma))**2)))

    @property
    def t(self):
        return ((self.depth/((1 + self.nu) * self.gamma**2))
               * (self.Edef['E1'] * (2*self.gamma * np.sinh(2*self.gamma) 
               - 4 * (self.gamma**2)) + (self.Edef['E2'] - self.Edef['E1']) 
               * (np.cosh(2*self.gamma) - 1 - 2*self.gamma**2))
               / (16 * ((np.sinh(self.gamma))**2)))

    def get_gamma(self, displacements, coords, k, t):

        return self.depth * np.sqrt(((1 - 2*self.nu) 
               * (np.trapz([theta*theta for theta in displacements['rotations']], coords)
               + 0.5 * np.sqrt(k / t) * ((displacements['vertical_displacements'][0])**2 
               + (displacements['vertical_displacements'][-1])**2))) / ((2 * (1 - self.nu))
               * (np.trapz([w*w for w in displacements['vertical_displacements']], coords)
               + 0.5 * np.sqrt(t / k) * ((displacements['vertical_displacements'][0])**2 
               + (displacements['vertical_displacements'][-1])**2))))

    def __str__(self):
        return f"Edef: {self.Edef} \nnu: {self.nu} \ndepth: {self.depth} \ngamma: {self.gamma} \nk: {self.k} \nt: {self.t}"
