import numpy as np
import compressorParameters as P
import pandas as pd
from collections import OrderedDict

class compressorStageData:
    parameters = [val.strip() for val in "beta_1, beta_2, alpha_1, alpha_2, C_w1, C_w2, C_1, C_2, V_w1, V_w2, V_1, V_2, U, M_1, M_2, Lambda, deHaller, r".split(",")]
    max_len = max([len(name) for name in parameters])
    
    def __init__(self, input_data, name) -> None:
        self.name = name
        self.data = OrderedDict()
        for param in self.parameters:
            if param in input_data:
                self.data[param] = input_data[param]
            else:
                self.data[param] = None
    
    def getFormattedColumns(self):
        names = {}
        specials = ["beta", "alpha", "Lambda", "deHaller"]
        for param in self.parameters:
            if "_" in param:
                temp = param[:param.find("_")+1] + "{" + param[param.find("_")+1:] + "}"
                names[param] = f"$\{temp}$" if param[:param.find("_")] in specials else f"${temp}$"
            else:
                names[param] = f"$\{param}$" if param in specials else f"${param}$"

        # names["Ca"] = "$C_a$"
        return names

    def toDataFrame(self):
        return pd.DataFrame(self.data, index=[self.name])

    def print(self):
        print(2*" " + self.name)
        for param in self.data:
            if self.data[param] is not None:
                print(4*" " + param + (self.max_len-len(param))*" " + ": {:.4f}".format(self.data[param]))
            else:
                print(4*" " + param + (self.max_len-len(param))*" " +": None")
        

class freeVortexCompressorMeanLine:
    gamma    = P.gamma
    R        = P.R
    cp       = P.cp
    Mach_max = P.Mach_max
    de_Haller_min = P.de_Haller_min

    beta_1m, beta_2m, alpha_1m, alpha_2m, Ca, Lambda_m, r_m, T_o1, N, delta_T_o, wdf = None, None, None, None, None, None, None, None, None, None, None
    
    def __init__(self, beta_1m, beta_2m, alpha_1m, alpha_2m, Ca, Lambda_m, r_m, T_o1, N, delta_T_o, wdf) -> None:
        for val in "beta_1m, beta_2m, alpha_1m, alpha_2m, Ca, Lambda_m, r_m, T_o1, N, delta_T_o, wdf".split(","):
            val = val.strip()
            exec(f"self.{val} = {val}")
        
        self.C_w1m = self.Ca*np.tan(self.alpha_1m)
        self.C_w2m = self.Ca*np.tan(self.alpha_2m)

    def calculate(self, r, name=None) -> compressorStageData:
        R = r/self.r_m
        U = 2*np.pi*r*self.N

        C_w1 = self.C_w1m * self.r_m / r
        # C_w2 = self.C_w2m * self.r_m / r

        # alpha_1 = np.rad2deg(np.arctan2(C_w1, self.Ca))
        # alpha_2 = np.rad2deg(np.arctan2(C_w2, self.Ca))
        
        # beta_1  = np.rad2deg(np.arctan2(V_w1, self.Ca))
        # beta_2  = np.rad2deg(np.arctan2(V_w2, self.Ca))
        
        Lambda = 1 - (1 - self.Lambda_m)/R**2
        # We already know lamda
        # Equations to help solve for betas
        eq1= (self.cp*self.delta_T_o) / (self.wdf*U*self.Ca) # = tan(B1) - tan(B2)
        eq2= 2*Lambda*U/self.Ca # = tan(B1) + tan(B2)
        
        # Solve for gas angles
        beta_1 = np.arctan((eq1+eq2)/2)
        beta_2 = np.arctan(eq2 - np.tan(beta_1))
        alpha_1 = np.arctan(U/self.Ca  - np.tan(beta_1))
        alpha_2 = np.arctan(U/self.Ca  - np.tan(beta_2))

        dCw = self.cp*self.delta_T_o / (self.wdf*U) # Cw2 - Cw1 = dCw
        C_w2 = dCw + C_w1
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


        # Solve for whirl velocities
        C_w1 = self.Ca*np.tan(alpha_1)
        C_w2 = self.Ca*np.tan(alpha_2)

        V_w1    = (U - C_w1)
        V_w2    = (U - C_w2)

        C_1  = np.sqrt(C_w1**2 + self.Ca**2)
        C_2  = np.sqrt(C_w2**2 + self.Ca**2)
        
        T_1  = self.T_o1 - (C_1**2)/(2*self.cp)
        T_2  = self.T_o1 - (C_2**2)/(2*self.cp)
    
        V_1  = np.sqrt((2*np.pi*self.N*r - C_w1)**2 + self.Ca**2)
        V_2  = np.sqrt((2*np.pi*self.N*r - C_w2)**2 + self.Ca**2)
        
        M_1  = V_1 / np.sqrt(self.gamma*self.R*T_1)
        M_2  = V_2 / np.sqrt(self.gamma*self.R*T_2)
        
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
        if Lambda_test < 0:
            if name is None:
                raise ValueError('degree of reaction check failed: Lambda = {:.4f}'.format(Lambda_test))
            else:
                raise ValueError('degree of reaction check failed at {}: Lambda = {:.4f}'.format(name, Lambda_test))
        
        # Check de Haller Criteria
        deHaller = V_2/V_1
        if deHaller < self.de_Haller_min:
            if name is None:
                raise ValueError('de Haller check failed: V2/V1 = {:.4f}'.format(deHaller))
            else:
                raise ValueError('de Haller check failed at {}: V2/V1 = {:.4f}'.format(name, deHaller))
        Ca = self.Ca

        beta_1 = np.rad2deg(beta_1)
        beta_2 = np.rad2deg(beta_2)
        alpha_1 = np.rad2deg(alpha_1)
        alpha_2 = np.rad2deg(alpha_2)

        return compressorStageData(locals(), name)