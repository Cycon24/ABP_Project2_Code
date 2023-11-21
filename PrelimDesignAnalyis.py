# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 00:17:05 2023

@author: cycon
"""

import EngineModule as EM
from parameters.OveralAnalysis import turbofan_kwargs 

# =============================================================================
#   Engine Input Parameters
# =============================================================================
Turbofan = EM.Turbofan_DoubleSpool(**turbofan_kwargs)
# Force full expansion
Turbofan.BP_nozzle.Pe = turbofan_kwargs['Pa']
Turbofan.nozzle.Pe = turbofan_kwargs['Pa']

# =============================================================================
#     Calculation
# =============================================================================
# Print Inputs
Turbofan.printInputs()
# Calculate 
Turbofan.calculate(True)