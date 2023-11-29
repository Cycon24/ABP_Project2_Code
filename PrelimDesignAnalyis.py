# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 00:17:05 2023

@author: cycon
"""
import numpy as np
import EngineModule as EM
from parameters.OveralAnalysis import turbofan_kwargs 
import warnings, pandas
from freeVortexCompressor import compressorStageData
warnings.filterwarnings("error")

# Adjustable parameters - Compressor
nc_p = 0.92 # >= 0.90 - Polytropic efficiency of overall compressor
Mt_1 = 1.19 # <= 1.20 - tip mach number rel to blade for first stage
Ca1  = 230  # m/s - Assuming Constant
Cw1  = 0    # m/s
N    = 200  # rev/s
dTos_ie = 45 # K - stag temp rise for first and last stages
num_extra_stages = 1 # The number of stages to add after initial #stages estimation
Lamda_after2nd = 0.5 # Degree of reaction at m line for stages 3+

# Adjustable parameters - Turbine
num_turb_stages = 2
phi = 0.78 # >= .78 limit for Flow Coefficient
psi = 3.3  # <= 3.3 limit for Blade Loading Coefficient
Lamda_turb = .50 # Chosen Degree of Reaction

# Gas Properties
R =   287 # J/kg*K
    # Air 
cp = 1005 # J/kg*K
gam = 1.4
    # Gas
cpg = 1148 # J/kg*K
gam_g = 4/3 

# Print Inputs for documentation
print('Controllable Inputs')
print('\tStage Isentropic Efficiency: ',nc_p)
print('\tMach number at tip of Stage1: ',Mt_1)
print('\tCompressor Inlet Axial Vel: ', Ca1, ' m/s')
print('\tCompressor Inlet Swirl Vel: ', Cw1, ' m/s')
print('\tCompressor rotational speed: ', N, ' rev/s')
print('\tStagnation temp rise of 1st and last stages: ', dTos_ie, ' K')
print('\tStages added onto approximation: ',num_extra_stages) 
print('\tAssumed Deg. of React. for stages 3+: ', Lamda_after2nd)
print('\n')
# =============================================================================
#     Calculation of Overall Engine Characteristics
# =============================================================================
turbofan_kwargs['npc'] = nc_p 
Turbofan = EM.Turbofan_DoubleSpool(**turbofan_kwargs)
# Force full expansion
Turbofan.BP_nozzle.Pe = turbofan_kwargs['Pa']
Turbofan.nozzle.Pe = turbofan_kwargs['Pa']

# Print Inputs
Turbofan.printInputs()
print('----------------------')
# Calculate 
print('\nOverall Enginge Analyis Outputs: ')
Turbofan.calculate(False)
Turbofan.fan.printOutputs()
Turbofan.HP_comp.printOutputs()
print('----------------------')


# =============================================================================
#   Compressor Stage Approximations
# =============================================================================
# Retrieve the compressor inlet and outlet conditions 
To1 = Turbofan.HP_comp.Toi
Po1 = Turbofan.HP_comp.Poi
Toe = Turbofan.HP_comp.Toe
Poe = Turbofan.HP_comp.Poe 
mdot_c = Turbofan.HP_comp.m_dot

dTo_comp = Toe - To1


# Calculate the mean radius and rotational velocity at mean line 
# based on a chosed tip mach number

# Static properties
T1   = To1 - (Ca1**2) / (2*cp)
P1   = Po1*(T1/To1)**(gam/(gam-1))
rho1 = P1/(R*T1)

# Find first radii
V1_t = Mt_1*np.sqrt(gam*R*T1)
U1_t = np.sqrt(V1_t**2 - Ca1**2)
r1_t = U1_t / (2*np.pi*N)
r1_r = r1_t *np.sqrt(1 - mdot_c/(rho1*Ca1*np.pi*r1_t**2))
r_m = (r1_t + r1_r)/2  # Constant through the stages

# Stage number estimation
U_m = 2*np.pi*r_m*N
beta_1m = np.arctan2(U_m, Ca1) # Radians
V_1m = Ca1/np.cos(beta_1m)
V_2m = 0.65*V_1m # Using de Haller crit
beta_2m = np.arccos(Ca1/(V_2m)) # Radians

# IGNORING WorkDoneFactor for initial estimations
dTo_s_estimate = U_m*Ca1*(np.tan(beta_1m)-np.tan(beta_2m)) / cp # K

stage_approx = dTo_comp / dTo_s_estimate
num_stages_approx = int(np.ceil(stage_approx))
print('Initial Aproximation for number of Stages:')
print('{:6.3f} => {} stages'.format(stage_approx, num_stages_approx))
print(f'Increasing to {num_stages_approx + num_extra_stages} stages')
num_stages_approx += num_extra_stages # add stagse to account for WDF losses in estimate

# Assuming a dTos = 30K for the first and last stages
dTo_s_est_stage = dTo_comp / num_stages_approx # All stages if even dTo distribution
dTo_s_est_remaining = (dTo_comp - 2*dTos_ie) / (num_stages_approx-2) # For all intermediate stages

dTos_int = dTo_s_est_remaining # for intermediate stages

# =============================================================================
# Compressor Stage Calculations
# =============================================================================
allStage_inputs = {
    'Ca': Ca1,
    'r_m':r_m,
    'nc_s': nc_p,
    'U_m': U_m,
    'mdot_c': mdot_c,
    'N': N
}

# Add all stages into an array to iterate through
stages = []
for i in range(0,num_stages_approx):
    if i == 0:
        # Need none here to calc own degree of reaction
        stage = EM.compressor_stage(i+1, None, dTos_ie, Cw1=Cw1, Toi=To1, Poi=Po1, **allStage_inputs)
    elif i == 1:
        # need 0 input as reminder to manually average DoR from stage 1
        stage = EM.compressor_stage(i+1, 0, dTos_int, **allStage_inputs)
    
    elif i == num_stages_approx - 1: # Last stage
        stage = EM.compressor_stage(i+1, Lamda_after2nd, dTos_ie, **allStage_inputs)
    
    else: # All intermediate stages
        # Assumes DoR = 0.5
        stage = EM.compressor_stage(i+1, Lamda_after2nd, dTos_int, **allStage_inputs)
    
    stages.append(stage)


stage_datas = []

# Propogate calculations through stages
for i, stage in enumerate(stages):
    # Check if in second stage and average lamda
    if i == 1:
        lam1 = stages[i-1].Lamda 
        lam3 = Lamda_after2nd 
        lam2 = (lam1 + lam3)/2
        stage.Lamda = lam2
    
    M_flag, deHal_flag = stage.calculate()
    if M_flag or deHal_flag:
        break
    if i != num_stages_approx -1:
        stage.forward(stages[i+1])
    # stage.printVelocityTrianges()
    stage_datas.append(stage.data)

for_formatting = compressorStageData({}, "")
with open("out.tex", 'w') as f:
    out =pandas.concat(stage_datas)
    names = for_formatting.getFormattedColumns()
    out = out.rename(columns=names)
    f.write(out.to_latex())

print('Finished calculations')   


# =============================================================================
# Turbine Calculations  
# =============================================================================
# Retrieve Tubine inlet and outlet values
Toi_turb = Turbofan.HP_turb.Toi
Toe_turb = Turbofan.HP_turb.Toe
Poi_turb = Turbofan.HP_turb.Poi
Poe_turb = Turbofan.HP_turb.Poe
dTo_turb = Toi_turb - Toe_turb

# Required Values
dTo_s_turb = dTo_turb/num_turb_stages # K
nt_p = turbofan_kwargs['npt']

# Thus, presure Ratio for per stage
Po4cycle = Poe_turb # 4 represents exit of stage 2, 3 is exit stage 1
Po3_over_Po4cycle = (1 - dTo_s_turb/Toi_turb)**(gam_g/(nt_p*(gam_g - 1)))
Po3 = Po4cycle * Po3_over_Po4cycle # Stage pressure between stage 1 and 2 of turb

# Velocities at Mean
Um_turb = np.sqrt(2*cpg*dTo_s_turb/psi)
Ca_turb = phi * Um_turb

# Gas angle and Swirl Angles Respectively
beta2  = np.arctan(0.5/phi * (psi/2 - 2*Lamda_turb)) # rad
beta3  = np.arctan(0.5/phi * (psi/2 + 2*Lamda_turb)) # rad
alpha2 = np.arctan(np.tan(beta2) + 1/phi) # rad
alpha3 = np.arctan(np.tan(beta3) - 1/phi) # rad

# Mean Radius
rm = Um_turb / (2*np.pi*N) # m

# Obtaining Bade Height
C3          = Ca_turb/np.cos(alpha3) # m/s
To3         = To_turb_in + dTos_per_stage # K
T3          = To3 - C3**2 / (2*cpg) # K
P3_over_Po3 = (T3/To3)**(gamma_g/(gamma_g - 1))
P3          = Po3 * P3_over_Po3 # bar
rho3        = P3*1e5 / (R*T3) # kg/m^3
h3          = mdot_a / (2*np.pi*rho3*Ca*rm) # m
rt3         = rm + h3/2 # m
rr3         = rm - h3/2 # m
Cw2m        = cpg * dTos_per_stage / Um # m/s
Cw2t        = Cw2m * rm/rt # m/s
C2t         = np.sqrt(Cw2t**2 + Ca**2) # m/s
Ut2         = Um * rt/rm # m/s
V2t         = np.sqrt((Cw2t - Ut2)**2 + Ca**2) # m/s
beta2t      = np.arccos(Ca / V2t) # rad
alpha2t     = np.arctan(Cw2t / Ca) # rad
T2t         = To1 - C2t**2 / (2*cpg) # K
M2t         = V2t / np.sqrt(gamma_g * R * T2t)\
