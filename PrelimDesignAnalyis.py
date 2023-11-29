# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 00:17:05 2023

@author: cycon
"""
import numpy as np
import EngineModule as EM
from parameters.OveralAnalysis import turbofan_kwargs 
import compressorParameters as comp
import warnings, pandas
from freeVortexCompressor import compressorStageData
import matplotlib.pyplot as plt
warnings.filterwarnings("error")

# Adjustable parameters
nc_p = 0.92 # >= 0.90 - Polytropic efficiency of overall compressor
Mt_1 = 1.19 # <= 1.20 - tip mach number rel to blade for first stage
Ca1  = 230 # m/s - Assuming Constant
Cw1  = 0    # m/s
N    = 200  # rev/s
dTos_ie = 45 # K - stag temp rise for first and last stages
num_extra_stages = 1 # The number of stages to add after initial #stages estimation
Lamda_after2nd = 0.5 # Degree of reaction at m line for stages 3+

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
cp = 1005 # J/kg*K
R =   287 # J/kg*K
gam = 1.4

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





out =pandas.concat(stage_datas)
for_formatting = compressorStageData({}, "")
with open("out.tex", 'w') as f:
    angles = out[["alpha_1", "alpha_2", "beta_1", "beta_2"]]
    absolute_velocities = out[["C_w1", "C_w2", "C_1", "C_2", "U"]]
    relative_velocities = out[["V_w1", "V_w2", "V_1", "V_2"]]
    misc = out[['M_1', 'M_2', 'Lambda', 'deHaller', 'r']]

    names = for_formatting.getFormattedColumns()
    angles = angles.rename(columns=names)
    absolute_velocities = absolute_velocities.rename(columns=names)
    relative_velocities = relative_velocities.rename(columns=names)
    misc = misc.rename(columns=names)
    
    f.write(
        "\\begin{center}\n" + angles.to_latex(float_format='%.3f') + "\\end{center}\n\n" +
        "\\begin{center}\n" + absolute_velocities.to_latex(float_format='%.3f') + "\\end{center}\n\n" +
        "\\begin{center}\n" + relative_velocities.to_latex(float_format='%.3f') + "\\end{center}\n\n" +
        "\\begin{center}\n" + misc.to_latex(float_format='%.3f')+ "\\end{center}\n\n" 
    )

out_T = out.T
# out.plot()
# plt.show()
# #print(out["beta_1"])
plt.rcParams["text.usetex"] = True

plt.figure()
for name in ["root", "mean", "tip"]:
    count = 1
    mach_data = []
    deHallar_data = []
    Lambda_data = []
    for item in out_T:
        if item[1] == name:
            mach_data.append(out_T[item[0]][name]["M_1"])
            count+=1
    plt.plot(list(range(1, count)), mach_data, label=name[0].upper() + name[1:])

plt.plot(list(range(1, count)), (count-1)*[comp.Mach_max], "--", label="Limit")
plt.title(f"Mach Number along the Stages")
plt.legend()
plt.ylabel("Mach number")
plt.xlabel("Stage")
plt.grid(True)
plt.minorticks_on()
plt.tight_layout()
plt.savefig(f"mach.png")

plt.figure()
for name in ["root", "mean", "tip"]:
    count = 1
    deHallar_data = []
    for item in out_T:
        if item[1] == name:
            deHallar_data.append(out_T[item[0]][name]["deHaller"])
            count+=1
            # print(out_T[item[0]][name]["M_1"])
    plt.plot(list(range(1, count)), deHallar_data, label=name[0].upper() + name[1:])

plt.plot(list(range(1, count)), (count-1)*[comp.de_Haller_min], "--", label="Limit")
plt.title(f"de-Haller Number along the Stages")
plt.legend()
plt.ylabel("de-Haller Number")
plt.xlabel("Stage")
plt.grid(True)
plt.minorticks_on()
plt.tight_layout()
plt.savefig(f"deHaller.png")

plt.figure()
for name in ["root", "mean", "tip"]:
    count = 1
    Lambda_data = []
    for item in out_T:
        if item[1] == name:
            Lambda_data.append(out_T[item[0]][name]["Lambda"])
            count+=1
    plt.plot(list(range(1, count)), Lambda_data, label=name[0].upper() + name[1:])

plt.title(f"Degree of Reaction $\Lambda$ along the Stages")
plt.legend()
plt.ylabel("Degree of Reaction $\Lambda$")
plt.xlabel("Stage")
plt.grid(True)
plt.minorticks_on()
plt.tight_layout()
plt.savefig(f"Lambda.png")
    # print(item)

print('Finished calculations')

