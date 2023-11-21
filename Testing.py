# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 14:23:02 2023

@author: cycon
"""

import EngineModule as EM


def Ex1():
    npc = 0.87 
    npt = 0.87 
    ni  = 0.95 
    nj  = 0.95 
    nm  = 0.99 
    dP_b = 0.06 
    nb  = 0.97 
    
    Ta = 242.7 # K
    Pa = 0.411 # bar
    rc = 8.0 
    Ttoi = 1200 # K - Turbine inlet temp
    
    mdot = 15 # kg/s
    Ci  = 260 # m/s
    
    gen_kwargs = {
        'Ta': Ta,
        'Pa': Pa,
        'Vinf':Ci,
        'rc':rc}
    
    intk = EM.Intake(**gen_kwargs, ni=ni,m_dot=mdot)
    comp = EM.Compressor(**gen_kwargs, np=npc)
    combust = EM.Combustor(**gen_kwargs, Toe=Ttoi, dPb_dec=dP_b, ni=nb)
    turb = EM.Turbine(comp, nm=nm, np=npt)
    exh  = EM.Nozzle(**gen_kwargs, ni=nj)
    
    intk.calculate()
    intk.forward(comp)
    comp.calculate()
    comp.forward(combust)
    combust.calculate()
    combust.forward(turb)
    turb.calculate()
    turb.forward(exh)
    exh.calculate()
    
    intk.printOutputs()
    comp.printOutputs()
    combust.printOutputs()
    turb.printOutputs()
    exh.printOutputs()
    
# =============================================================================
#     Exammple 2
# =============================================================================
np = 0.90
ni  = 0.95 

dP_b = 0.06 

Ta = 216.66 # K
Pa = 0.224 # bar
rfan = 1.55
rc = 35/rfan
Ttoi = 1350 # K - Turbine inlet temp
BPR = 6.2

mdot = 220 # kg/s
Mi  = 0.85 # Mach

gen_kwargs = {
    'Ta': Ta,
    'Pa': Pa,
    'Minf':Mi,
    'np':np}

intk = EM.Intake(**gen_kwargs, ni=ni,m_dot=mdot)
fan = EM.Compressor(**gen_kwargs, rc=rfan, BPR=BPR)
hpc = EM.Compressor(**gen_kwargs, rc=rc)
combust = EM.Combustor(**gen_kwargs, Toe=Ttoi, dPb_dec=dP_b)
hpt = EM.Turbine(hpc, **gen_kwargs)
lpt = EM.Turbine(fan, **gen_kwargs)
exh_c  = EM.Nozzle('cold', **gen_kwargs)
exh_h  = EM.Nozzle(**gen_kwargs)

intk.calculate()
intk.forward(fan)
fan.calculate()
fan.forward(hpc,exh_c)
# Cold
exh_c.calculate()
# Hot
hpc.calculate()
hpc.forward(combust)
combust.calculate()
combust.forward(hpt)
hpt.calculate()
hpt.forward(lpt)
lpt.calculate()
lpt.forward(exh_h)
exh_h.calculate()

intk.printOutputs()
fan.printOutputs()
exh_c.printOutputs()
hpc.printOutputs()
combust.printOutputs()
hpt.printOutputs()
lpt.printOutputs()
exh_h.printOutputs()