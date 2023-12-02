import numpy as np
import warnings
warnings.filterwarnings("error") # switch warnings to error, catches floating point errors

from Compressor import ConstantAxialVelocityConstantMeanLineCompressor      as Compressor
from Compressor import DummyConstantAxialVelocityConstantMeanLineCompressor as DummyCompressor
from Compressor import drawCompressorGeometry                               as drawCompressor
from Compressor import outputData                                           as outputCompressorData

from Stator     import ConstantAxialVelocityConstantMeanLineStator  as Stator
from Stator     import outputData                                   as outputStatorData

from Turbine    import ConstantAxialVelocityConstantMeanLineTurbine as Turbine
from Turbine    import drawTurbineGeometry                          as drawTurbine
from Turbine    import outputData                                   as outputTurbineData

from tools      import drawTotalGeometry


gamma   = 1.4         # gamma for air
c_pa    = 1.005*10**3 # cp for air
c_pg    = 1.148*10**3 # cp for gas
R       = 287         # R for air and gas
Delta_H = 43100*10**3 # fuel heating value
eta_b   = 1           # assume perfect combustion

T_o4    = 2275.0 # combustion exit temperature
T_a     = 288.0  # Required ambient temperature
P_a     = 1.01e5 # Required ambient pressure (Pa)
n_inf_f = 0.92   # Fan Polytropic efficiency
pi_f    = 1.4    # Required Fan pressure ratio
pi_c    = 15     # Required Compressor pressure ratio
BPR     = 6      # Required Bypass Ratio
m_dot_a = 280    # Mass flow into the engine in kg/s


####
#  Adjustable parameters, we iterated on these
####
n_inf_c = 0.92  # >= 0.90 - Polytropic efficiency of overall compressor
M_1t    = 1.19  # <= 1.20 - tip mach number rel to blade for first stage, if the flow is axial
C_a1    = 231   # inlet axial velocity, m/s - Assuming Constant
C_w1t   = 50    # inlet swirl at the tip, m/s, reduces inlet Mach number
N       = 200   # rotation rate of rotor/shaft, rev/s

Delta_T_os_ie    = 45  # K - stag temp rise for first and last stages
num_extra_stages = 1   # The number of stages to add after initial stage count estimation
Lamda_after2nd   = 0.5 # Degree of reaction at mean line for stages 3+

# turbine, choose at the mean line
psi                 = 2.8  # blade loading coefficient
phi                 = 0.84 # flow coefficient
Lambda_t            = 0.5  # degree of reaction
num_turbine_stages  = 2   # number of stages that we desired


##########
#
# Compressor Calculations
#
##########


"""
Computes the mass flow through the compressor

"""
def getMassFlowAtCompressorInlet():
    return m_dot_a/(BPR + 1)

"""
Computes the stagnation properties (temperature and pressure)
at the compressor inlet

"""
def getStagnationPropertiesAtCompressorInlet():
    T_o2 = T_a # isentropic inlet and the engine is not moving
    P_o2 = P_a # isentropic inlet and the engine is not moving
    
    P_o13 = pi_f*P_o2 # stagnation pressure after the fan
    T_o13 = T_o2*pi_f**((gamma-1)/(gamma*n_inf_f)) # stagnation temp after the fan
    return T_o13, P_o13


"""
Computes the stagnation properties (temperature and pressure)
at the compressor exit

"""
def getStagnationPropertiesAtCompressorExit(T_inlet, P_inlet):
    P_o_exit = P_inlet*pi_c # stagnation pressure after the compressor

    # stagnation temp after the compressor
    T_o_exit = T_inlet*(P_o_exit/P_inlet)**((gamma - 1)/(gamma*n_inf_c))

    return T_o_exit, P_o_exit

"""
Computes the static properties (temperature and pressure)
at the compressor inlet, first computes inlet velocity, then finds
the static properties using the stagnation properties

"""
def getStaticPropertiesAtCompressorInletTip():
    # finding velocity
    C_1t = np.sqrt(C_w1t**2 + C_a1**2)

    # just applying stagnation equations and ideal gas law
    T_1t   = T_o1 - (C_1t**2) / (2*c_pa)
    P_1t   = P_o1*(T_1t/T_o1)**(gamma/(gamma-1))
    rho_1t = P_1t/(R*T_1t)

    return T_1t, P_1t, rho_1t

"""
Estimates the number of stages using primarily the chosen tip mach number
but the rotation rate also affects it

"""
def estimateNumberOfCompressorStages():
    V_1t = M_1t*np.sqrt(gamma*R*T_1t) # relative velocity from chosen Mach number
    V_w1t = np.sqrt(V_1t**2 - C_a1**2) # # tip relative swirl, from velocity triangle
    U_t  = V_w1t + C_w1t # tip tangential velocity, from velocity triangle
    r_t  = U_t / (2*np.pi*N) # tip radius from angular velocity equation
    r_r  = r_t *np.sqrt(1 - m_dot_c/(rho_1t*C_a1*np.pi*r_t**2)) # root radius from continuity
    r_m  = (r_t + r_r)/2  # mean radius, Constant through the stages

    # Stage number estimation
    C_w1m   = C_w1t*r_t/r_m # apply free vortex condition
    U_m     = 2*np.pi*r_m*N # mean tangential velocity, from velocity triangle
    beta_1m = np.arctan2(V_w1t, C_a1) # inlet gas angle at mean line, Radians
    V_1m    = C_a1/np.cos(beta_1m) # inlet relative velocity at mean line
    V_2m    = 0.65*V_1m # max deHaller criteria, to get most performance
    beta_2m = np.arccos(C_a1/(V_2m)) # exit gas angle at mean line, Radians

    # IGNORING WorkDoneFactor for initial estimations 
    # per stage change in stagnation temperature
    Delta_T_os_estimate = U_m*C_a1*(np.tan(beta_1m)-np.tan(beta_2m)) / c_pa # K

    # finding the change in stagnation temp across the compressor
    T_o_exit, _ = getStagnationPropertiesAtCompressorExit(T_o1, P_o1)
    Delta_T_o_comp = T_o_exit - T_o1

    # finding the number of stages need to reach the total change in stag. temp.
    stage_approx = Delta_T_o_comp / Delta_T_os_estimate
    num_stages_approx = int(np.ceil(stage_approx))
    
    return num_stages_approx, r_m, C_w1m

"""
computes the stag. temp. change in intermediate stages.
All stages, but first and last

"""
def getIntermediateStageTemperatureRise(num_stages):
    # finding the change in stagnation temp across the compressor
    T_o_exit, _ = getStagnationPropertiesAtCompressorExit(T_o1, P_o1)
    Delta_T_o_comp = T_o_exit - T_o1

    # finding the stag. temp change for the intermediate stages,
    # by first reducing by stag temp change in first and last stage
    Delta_T_os_intermediate = (Delta_T_o_comp - 2*Delta_T_os_ie) / (num_stages-2)

    return Delta_T_os_intermediate


# getting the mass flow in the compressor
m_dot_c            = getMassFlowAtCompressorInlet()

# getting the stag props at the inlet to the compressor
T_o1, P_o1         = getStagnationPropertiesAtCompressorInlet()

# getting the compressor exit stag. props.
T_o_exit, P_o_exit = getStagnationPropertiesAtCompressorExit(T_o1, P_o1)

# getting the static props at the inlet to the compressor
T_1t, P_1t, rho_1t    = getStaticPropertiesAtCompressorInletTip()

# getting the approximate number of stages required (no work done factor)
# also getting the mean line radius for use all the time
num_stages_approx, r_m, C_w1m = estimateNumberOfCompressorStages()

print('Initial Aproximation for number of Stages:')
print(f'  {num_stages_approx} stages')
print(f'  Increasing to {num_stages_approx + num_extra_stages} stages')

# adding the additional amount of stages chosen by us (to make sure the requirements
# are satisfied)
num_stages = num_stages_approx + num_extra_stages # add stagse to account for WDF losses in estimate

# getting the required intermediate stage stag. temp. rise
Delta_T_os_intermediate = getIntermediateStageTemperatureRise(num_stages)

# holds the info
compressors = []
compressor_data = []

# used in the loop
T_o_next = T_o1
P_o_next = P_o1

# iterating over the compressors
print("Compressors")
for i in range(num_stages):
    # handles the first stage compressor
    if i == 0:
        # setting the stag. temp change to the value for the first and last stage
        Delta_T_o_required = Delta_T_os_ie

        # this tells the code to find the degree of reaction
        Lambda = None

        # specifying the swirl
        _C_w1m = C_w1m
    
    # handles the second stage compressor
    elif i == 1:
        # average last stage Lambda and next stage Lambda to get the second
        # stage's Lambda, "bridges the gap"
        Lambda = (compressors[-1].mean["Lambda"] + Lamda_after2nd)/2

        # setting the stag. temp change to the value for the intermediate stages
        Delta_T_o_required = Delta_T_os_intermediate

        # setting this to None tells code to find it
        _C_w1m = None
    
    # handles the last stage compressor
    elif i == num_stages - 1:
        # sets degree of reaction to what we want after the second stage
        # this is 0.5
        Lambda = Lamda_after2nd

        # finding Delta_T required by using the pressure ratio and the polytropic
        # efficiency as the single stage efficiency
        Delta_T_o_required = T_o_next/n_inf_c*((P_o_exit/P_o_next)**((gamma-1)/gamma) - 1)

        # setting this to None tells code to find it
        _C_w1m = None
    
    # handles all other stages of the compressor
    else:
        # sets degree of reaction to what we want after the second stage
        # this is 0.5
        Lambda = Lamda_after2nd

        # setting the stag. temp change to the value for the intermediate stages
        Delta_T_o_required = Delta_T_os_intermediate

        # setting this to None tells code to find it
        _C_w1m = None

    # calling the compressor class, adds information for later calculation
    compressor = Compressor(
        number      = i+1,
        C_am        = C_a1,
        r_m         = r_m,
        n_inf_c     = n_inf_c,
        m_dot       = m_dot_c,
        Delta_T_o   = Delta_T_o_required,
        N           = N
    )

    # calculates the stage parameters at root mean and tip
    # also sets the next stage's inlet pressure and temperature
    T_o_next, P_o_next = compressor.calculate(T_o_next, P_o_next, Lambda, C_w1m=_C_w1m)

    # collecting the data
    compressors.append(compressor)
    compressor_data.append(compressor.getData())
print()

# dumps all of the data, not important in the flow of calculations
outputCompressorData(compressor_data, compressors)

# adding a dummy compressor to represent what flow conditions are required
# in the combustor, used for the last stage stator
combustor_inlet = DummyCompressor(C_a1, r_m, m_dot_c, T_o_next, P_o_next)
compressors.append(combustor_inlet)
compressor_data.append(combustor_inlet.getData())

# drawing the compressor geometry
drawCompressor(compressor_data)


##########
#
# Stator Calculations
#
##########


# holds info
stators = []
stator_data = []

# iterating over the stators
print("Stators")
for i in range(0, num_stages):
    # creates a stator that statisfies the previous stage's and the next stage's require
    # flow angles and speeds
    stator = Stator(i+1, compressors[i], compressors[i+1])

    # collecting the data
    stators.append(stator)
    stator_data.append(stator.getData())
print()

# outputs the stator data
outputStatorData(stator_data)


##########
#
# Turbine Calculations
#
##########


"""
calculates the fuel to air ratio, f, from inlet and outlet
temperatures of the combustor

"""
def calc_f(T_o_out, T_o_in):
    return (c_pg*T_o_out - c_pa*T_o_in)/(eta_b*(Delta_H - c_pg*T_o_out))


# getting fuel mass flow from fuel to air ratio, adding this to
# compressor mass flow to get turbine mass flow
m_dot_f = m_dot_c*calc_f(T_o4, T_o_exit)
m_dot_t = m_dot_f + m_dot_c

# get compressor change in stag. temp.
Delta_T_oc = T_o_exit - T_o1

# get the total change in stag. temp. for the turbine from a
# power balance
Delta_T_ot = (c_pa*Delta_T_oc*m_dot_c)/(c_pg*m_dot_t)

# get exit temperature
T_o5 = T_o4 - Delta_T_ot

# assume equal power out of each stage
Delta_T_os = Delta_T_ot/num_turbine_stages

# finding tangential velocity from blade loading coefficient
# at the mean line to later get C_a
U_m = np.sqrt(2*c_pg*Delta_T_os/psi)

# get C_a from flow coefficient, this is constant
C_a = phi*U_m

# finding the mean radius, this is constant
r_m = U_m/(2*np.pi*N)

print("Turbine")

# hold turbine info
turbines = []
turbine_data = []

# setting the inlet stag. properties to those at the exit of 
# the combustor
T_o_next = T_o4
P_o_next = P_o_exit # no pressure loss in combustor

# iterating over the stages
for i in range(num_turbine_stages):
    # calling the turbine class, the input information is used for later calculation
    stage = Turbine(
        number=i+1, C_am=C_a, r_m=r_m, m_dot=m_dot_t, Delta_T_o=Delta_T_os, Lambda_m=Lambda_t, phi_m=phi, psi_m=psi, N=N
    )

    # setting the stag. props. of the next stage to ones at the exit
    # of this stage
    T_o_next, P_o_next = stage.calculate(T_o_next, P_o_next)

    # remembering the data
    turbines.append(stage)
    turbine_data.append(stage.getData())
print()

# outputting the data to files and plots
outputTurbineData(turbine_data, turbines)

# drawing the turbines
drawTurbine(turbine_data)

# drawing the turbine and compressor together
drawTotalGeometry(compressor_data, turbine_data)