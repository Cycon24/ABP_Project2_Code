# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 20:08:05 2023

@author: cycon
"""
# Double Spool Turbofan
turbofan_kwargs = {
    'Ta': 288.0, # Req - Pulled from hand calculations
    'Pa': 1.01e5, # Req
    'Vinf': 0, # Or Minf
    'Minf': None, # Or Vinf, if none its assumed stationary
    'ni': 1, # Inlet
    'nj': 1, # Nozzle
    'nf': None, # Compressor - Isentropic
    'nc': None, # Compressor - Isentropic
    'nt': None, # Turbine - Isentropic
    'nt_lp': None, # LP Turbine - Isentropic
    'nb': 1.0, # Combustor
    'nm': 0.99, # Mechanical
    'npf': 0.92, # Fan - Polytropic
    'npc': 0.90, # Compressor - Polytropic
    'npt': 0.92, # Turbine - Polytropic
    'npt_lp': 1, # LP Turbine - Polytropic
    'dP_combustor': 0, # Decimal pressure drop in combustor
    'rfan': 1.4, # Req - Fan PR
    'rc': 15,   # Req - Compressor PR
    'BPR': 6,     # Req Bypass Ratio
    'T_turb_in': 2275, # Req - K - Turbine inlet temp
    'mdot_a': 280,   # kg/s
    'Q_fuel': 43100 # kJ/kg
    }