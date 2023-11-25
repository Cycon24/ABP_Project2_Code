# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 00:17:05 2023

@author: cycon
"""

import EngineModule as EM
from parameters.OveralAnalysis import turbofan_kwargs 

# =============================================================================
#     Calculation of Overall Engine Characteristics
# =============================================================================
Turbofan = EM.Turbofan_DoubleSpool(**turbofan_kwargs)
# Force full expansion
Turbofan.BP_nozzle.Pe = turbofan_kwargs['Pa']
Turbofan.nozzle.Pe = turbofan_kwargs['Pa']

# Print Inputs
Turbofan.printInputs()
# Calculate 
Turbofan.calculate(True)

# =============================================================================
#   Compressor Stage Calculations
# =============================================================================
# Take out the compressor inlet and outlet conditions 
To1 = Turbofan.HP_comp.Toi
Po1 = Turbofan.HP_comp.Poi
Toe = Turbofan.HP_comp.Toe
Poe = Turbofan.HP_comp.Poe 

