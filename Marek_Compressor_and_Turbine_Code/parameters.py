
'''
This class contains parameters for the compressor.

These parameters are stored in a class so that they can be passed
to other classes via inheritance, it makes the code very clean

'''
class CompressorBaseParameters:
    def __init__(self) -> None:
        self.gamma = 1.4
        self.c_p = 1.005*10**3
        self.R  = 287

        self.T_a = 288.0 # Req - Pulled from hand calculations
        self.P_a = 1.01e5 # Req
        self.n_inf_f = 0.92 # Fan - Polytropic
        self.pi_f = 1.4 # Req - Fan PR
        self.pi_c = 15   # Req - Compressor PR
        self.BPR = 6     # Req Bypass Ratio
        self.m_dot_a = 280   # kg/s

        self.de_Haller_min = 0.65
        self.Mach_max = 1.2


'''
This class contains parameters for the turbine.

These parameters are stored in a class so that they can be passed
to other classes via inheritance, it makes the code very clean

'''
class TurbineBaseParameters:
    def __init__(self) -> None:
        self.gamma  = 1.333  # gamma for gas
        self.c_p    = 1.148*10**3 # c_p for gas
        self.R      = 287

        self.n_inf_t = 0.92 # polytropic efficiency for the turbine

        self.phi_min = 0.78 # min flow coefficient
        self.psi_max = 3.3  # mas blade loading