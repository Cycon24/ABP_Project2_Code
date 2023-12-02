import numpy as np
from parameters  import TurbineBaseParameters

from dataClasses import TurbineStageData


"""
Class for a turbine stage with a free vortex design and a constant
mean line

"""
class freeVortexCompressorMeanLine(TurbineBaseParameters):
    # initializing vars
    Ca, Lambda_m, psi_m, phi_m, r_m, T_o1, N, Delta_T_o, P_o1 = None, None, None, None, None, None, None, None, None

    def __init__(self, Ca, Lambda_m, psi_m, phi_m, r_m, T_o1, N, Delta_T_o, P_o1) -> None:
        super().__init__()
        # setting all the class attributes
        for val in "Ca, Lambda_m, psi_m, phi_m, r_m, T_o1, N, Delta_T_o, P_o1".split(","):
            val = val.strip()
            exec(f"self.{val} = {val}")
    

        # mean line tangential velocity
        self.U_m = 2*np.pi*self.r_m*self.N

        # mean line gas angles from flow parameters
        self.beta_2m  = np.arctan(1/(2*self.phi_m)*(self.psi_m/2 - 2*self.Lambda_m))
        self.beta_3m  = np.arctan(1/(2*self.phi_m)*(self.psi_m/2 + 2*self.Lambda_m))

        # mean line absolute angles
        self.alpha_2m = np.arctan(np.tan(self.beta_2m) + 1/self.phi_m) # rad
        self.alpha_3m = np.arctan(np.tan(self.beta_3m) - 1/self.phi_m) # rad

        # mean line swirl velocities
        self.C_w2m    = self.Ca*np.tan(self.alpha_2m)
        self.C_w3m    = self.Ca*np.tan(self.alpha_3m)

    def calculate(self, r, name):

        # getting the absolute angles from free vortex conditions
        alpha_2 = np.arctan(self.r_m/r*np.tan(self.alpha_2m))
        alpha_3 = np.arctan(self.r_m/r*np.tan(self.alpha_3m))

        # getting the gas angles from free vortex conditions
        beta_2  = np.arctan(self.r_m/r*np.tan(self.alpha_2m) - (r/self.r_m)*self.U_m/self.Ca)
        beta_3  = np.arctan(self.r_m/r*np.tan(self.alpha_3m) + (r/self.r_m)*self.U_m/self.Ca)

        # inlet and outlet swirl from velocity triangles
        C_w2   = self.Ca*np.tan(alpha_2)
        C_w3   = self.Ca*np.tan(alpha_3)

        # inlet and outlet absolute velocity from velocity triangles
        C_2     = np.sqrt(self.Ca**2 + C_w2**2)
        C_3     = np.sqrt(self.Ca**2 + C_w3**2)

        # tangential velocity
        U = 2*np.pi*r*self.N

        # inlet and outlet relative velocity from velocity triangles
        V_2 = self.Ca/np.cos(beta_2)
        V_3 = self.Ca/np.cos(beta_3)

        # inlet and outlet relative swirl velocity from velocity triangles
        V_w2 = self.Ca*np.tan(beta_2)
        V_w3 = self.Ca*np.tan(beta_3)

        # flow properties from their definitions
        phi = self.Ca/U
        psi = 2*self.c_p*self.Delta_T_o/U**2
        Lambda = self.Ca/(2*U)*(np.tan(beta_3) - np.tan(beta_2))

        ####
        # performing checks
        ####
        if phi < self.phi_min:
            raise ValueError('Flow Coefficient check failed at {}: phi = {:.4f}'.format(name, phi))
        
        if psi > self.psi_max:
            raise ValueError('Blade loading check failed at {}: psi = {:.4f}'.format(name, psi))
    
        if Lambda < 0 or Lambda >= 1:
            raise ValueError('Degree of reaction check failed at {}: Lambda = {:.4f}'.format(name, Lambda))

        beta_2 = np.rad2deg(beta_2)
        beta_3 = np.rad2deg(beta_3)
        alpha_2 = np.rad2deg(alpha_2)
        alpha_3 = np.rad2deg(alpha_3)


        # returning the data in the form of a class
        return TurbineStageData(locals(), name) # saves all the data to a class and returns it

        