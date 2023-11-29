# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 17:14:10 2023

@author: cycon
"""
import numpy as np
import EngineModule as EM
from parameters.OveralAnalysis import turbofan_kwargs 
import compressorParameters as comp
import warnings, pandas
from freeVortexCompressor import compressorStageData
import matplotlib.pyplot as plt
# TURBINE TEST FILE

# Adjustable parameters
nc_p = 0.92 # >= 0.90 - Polytropic efficiency of overall compressor
Mt_1 = 1.19 # <= 1.20 - tip mach number rel to blade for first stage
Ca1  = 230 # m/s - Assuming Constant
Cw1  = 0    # m/s
N    = 200  # rev/s
dTos_ie = 45 # K - stag temp rise for first and last stages
num_extra_stages = 1 # The number of stages to add after initial #stages estimation
Lamda_after2nd = 0.5 # Degree of reaction at m line for stages 3+

# Turbine
psi_m1 = 3.3
phi_m1 = 0.78 
Lamda_turb = 0.5

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
#  Turbine Testing
# =============================================================================
turb = EM.turbine_free_vortex(Turbofan.HP_turb, 2, psi_m1, phi_m1, Lamda_turb, turbofan_kwargs['npt'], N)
print(turb.calculate())