# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 13:34:16 2023

@author: cycon
"""
from EngineModule import Turbofan_SingleSpool as Turbofan
import sys
sys.path.append('/parameters/')
import parameters.realCycle as RC 
import parameters.idealCycle as IC 
# =============================================================================
#   Engine Input Parameters
# =============================================================================

turbofan_kwargs = IC.turbofan_kwargs

# Print Inputs
print('Inputs')
for key,val in turbofan_kwargs.items():
    print('\t {}  =  {}'.format(key,val))
    
# =============================================================================
#     Calculation
# =============================================================================
Engine = Turbofan(**turbofan_kwargs)
# Force full expansion
Engine.BP_nozzle.Pe = turbofan_kwargs['Pa']
Engine.nozzle.Pe = turbofan_kwargs['Pa']
# Calculate 
Engine.calculate()