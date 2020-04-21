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
        return ((1-self.nu)/((1+self.nu)*(1-2*self.nu)*self.depth))*(self.Edef['E1']*(2*self.gamma*np.sinh(2*self.gamma)+4*(self.gamma**2))+(self.Edef['E2']-self.Edef['E1'])*(np.cosh(2*self.gamma)-1+2*self.gamma**2))/(8*((np.sinh(self.gamma))**2))

    @property
    def t(self):
        return (self.depth/((1+self.nu)*self.gamma**2))*(self.Edef['E1']*(2*self.gamma*np.sinh(2*self.gamma)-4*(self.gamma**2))+(self.Edef['E2']-self.Edef['E1'])*(np.cosh(2*self.gamma)-1-2*self.gamma**2))/(16*((np.sinh(self.gamma))**2))


    def get_gamma(self, displacements, coords):
        return self.depth*np.sqrt(((1-2*self.nu)/(2*(1-self.nu)))*(np.trapz([theta*theta for theta in displacements['rotations']], coords)/np.trapz([w*w for w in displacements['vertical_displacements']], coords)))
#         function y=GAMMA(Translation,Rotation,Coords,nu,H)

# y=H*sqrt(((1-2*nu)/(2*(1-nu)))*(trapz(Rotation.^2,Coords)/trapz(Translation.^2,Coords))); [w*w for w in vertical_displacements]
    
    def __str__(self):
            return f"Edef: {self.Edef} \nnu: {self.nu} \ndepth: {self.depth} \ngamma: {self.gamma} \nk: {self.k} \nt: {self.t}"

# function y=VlasovParam(Foundation, Gamma)

# E=Foundation.Young;
# nu=Foundation.Poisson;
# H=Foundation.Depth;
# E1=E(1);
# E2=E(2);
# k=((1-nu)/((1+nu)*(1-2*nu)*H))*...
#     (E1*(2*Gamma*np.sinh(2*Gamma)+4*(Gamma**2))+(E2-E1)*(np.cosh(2*Gamma)-1+2*Gamma**2))/(8*((sinh(Gamma))^2));


# t=(H/((1+nu)*Gamma^2))*...
#     (E1*(2*Gamma*sinh(2*Gamma)-4*(Gamma^2))+(E2-E1)*(np.cosh(2*Gamma)-1-2*Gamma^2))/(16*((sinh(Gamma))^2));

# y=[k ; t];

