import numpy as np
from parameters  import CompressorBaseParameters
from dataClasses import CompressorStageData

"""
Class for a compressor stage with a free vortex design and a constant
mean line

"""
class freeVortexCompressorMeanLine(CompressorBaseParameters):
    # initializing vars
    beta_1m, beta_2m, alpha_1m, alpha_2m, Ca, Lambda_m, r_m, T_o1, N, Delta_T_o, wdf, n_inf_c, P_o1 = None, None, None, None, None, None, None, None, None, None, None, None, None

    def __init__(self, beta_1m, beta_2m, alpha_1m, alpha_2m, Ca, Lambda_m, r_m, T_o1, P_o1, N, Delta_T_o, wdf, n_inf_c) -> None:
        super().__init__()

        # setting all the class attributes
        for val in "beta_1m, beta_2m, alpha_1m, alpha_2m, Ca, Lambda_m, r_m, T_o1, N, Delta_T_o, wdf, n_inf_c, P_o1".split(","):
            val = val.strip()
            exec(f"self.{val} = {val}")
        
        # getting swirl velocities at inlet and outlet
        self.C_w1m = self.Ca*np.tan(self.alpha_1m)
        self.C_w2m = self.Ca*np.tan(self.alpha_2m)

    def calculate(self, r, name=None) -> CompressorStageData:
        # find stag props
        P_o3 = self.P_o1*(1 + self.n_inf_c*self.Delta_T_o/self.T_o1)**(self.gamma/(self.gamma - 1))
        T_o3 = self.T_o1 + self.Delta_T_o
    
        # get tangential velocity
        U = 2*np.pi*r*self.N

        # getting inlet swirl velocity        
        C_w1 = self.C_w1m * self.r_m / r

        # getting degree of reaction
        Lambda = 1 - (1 - self.Lambda_m)/(r/self.r_m)**2

        # Equations to help solve for betas
        eq1= (self.c_p*self.Delta_T_o) / (self.wdf*U*self.Ca) # = tan(B1) - tan(B2)
        eq2= 2*Lambda*U/self.Ca # = tan(B1) + tan(B2)
        
        # Solve for gas angles
        beta_1 = np.arctan((eq1+eq2)/2)
        beta_2 = np.arctan( eq2 - np.tan(beta_1))

        # get absolute angles
        alpha_1 = np.arctan(U/self.Ca  - np.tan(beta_1))
        alpha_2 = np.arctan(U/self.Ca  - np.tan(beta_2))

        # getting rotor exit swirl
        dCw = self.c_p*self.Delta_T_o / (self.wdf*U) # Cw2 - Cw1 = dCw
        C_w2 = dCw + C_w1

        # Solve for whirl velocities
        C_w1 = self.Ca*np.tan(alpha_1)
        C_w2 = self.Ca*np.tan(alpha_2)

        # relative swirl velocities
        V_w1    = (U - C_w1)
        V_w2    = (U - C_w2)

        # absolute velocities
        C_1  = np.sqrt(C_w1**2 + self.Ca**2)
        C_2  = np.sqrt(C_w2**2 + self.Ca**2)
        
        # static temperatures
        T_1  = self.T_o1 - (C_1**2)/(2*self.c_p)
        T_2  = T_o3 - (C_2**2)/(2*self.c_p)
    
        # get the relative velocity magnitudes (from velocity triangles)
        V_1  = np.sqrt(V_w1**2 + self.Ca**2)
        V_2  = np.sqrt(V_w2**2 + self.Ca**2)
        
        # getting relative mach numbers
        M_1  = V_1 / np.sqrt(self.gamma*self.R*T_1)
        M_2  = V_2 / np.sqrt(self.gamma*self.R*T_2)
        
        ###
        # performing checks, on the design parameters
        ###
        if M_1 > self.Mach_max:
            if name is None:
                raise ValueError('Max Mach exceeded: M_1 = {:.4f}'.format(M_1))
            else:
                raise ValueError('Max Mach exceeded at {} : M_1 = {:.4f}'.format(name, M_1))

        if M_2 > self.Mach_max:
            if name is None:
                raise ValueError('Max Mach exceeded: M_2 = {:.4f}'.format(M_2))
            else:
                raise ValueError('Max Mach exceeded at {} : M_2 = {:.4f}'.format(name, M_2))

        Lambda_test = self.Ca/(2*U)*(np.tan(beta_1) + np.tan(beta_2))
        if Lambda_test < 0 or Lambda_test >= 1:
            if name is None:
                raise ValueError('degree of reaction check failed: Lambda = {:.4f}'.format(Lambda_test))
            else:
                raise ValueError('degree of reaction check failed at {}: Lambda = {:.4f}'.format(name, Lambda_test))
        
        if abs(Lambda - Lambda_test) > 1e-14:
            raise ValueError(f"Lambda did not match its alternative calculation: given = {Lambda},  alt = {Lambda_test}")

        # Check de Haller Criteria
        deHaller = V_2/V_1
        if deHaller < self.de_Haller_min:
            if name is None:
                raise ValueError('de Haller check failed: V2/V1 = {:.4f}'.format(deHaller))
            else:
                raise ValueError('de Haller check failed at {}: V2/V1 = {:.4f}'.format(name, deHaller))

        # converting angles to degrees, for printing later
        beta_1 = np.rad2deg(beta_1)
        beta_2 = np.rad2deg(beta_2)
        alpha_1 = np.rad2deg(alpha_1)
        alpha_2 = np.rad2deg(alpha_2)

        # returning the data in the form of a class
        return CompressorStageData(locals(), name) # saves all the data to a class and returns it
    

"""
Depricated:
# Solve for Gas angles
# beta_1_test = np.arctan((U - C_w1) / self.Ca)
# beta_2_test = np.arctan((U - C_w2) / self.Ca)
# alpha_1_test = np.arctan(C_w1/self.Ca) # Should be 0 for first stage
# alpha_2_test = np.arctan(C_w2/self.Ca)

# if alpha_1 != alpha_1_test:
#     raise ValueError("alpha_1s not equal")

# if alpha_2 != alpha_2_test:
#     raise ValueError("alpha_2s not equal")

# if beta_1 != beta_1_test:
#     raise ValueError("beta_1s not equal")

# if beta_2 != beta_2_test:
#     raise ValueError(f"beta_2s not equal: {beta_2}, {beta_2_test}")

# C_w2 = self.C_w2m * self.r_m / r

# alpha_1 = np.rad2deg(np.arctan2(C_w1, self.Ca))
# alpha_2 = np.rad2deg(np.arctan2(C_w2, self.Ca))

# beta_1  = np.rad2deg(np.arctan2(V_w1, self.Ca))
# beta_2  = np.rad2deg(np.arctan2(V_w2, self.Ca))
"""