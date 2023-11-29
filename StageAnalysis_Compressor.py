# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 14:52:26 2023

@author: cycon
"""
import numpy as np
from WorkDoneFactor import interp_wdf

# Inputs / Knobs to turn
nc_p = 0.92 # >= 0.90
Mt   = 1.15 # <= 1.20
Ca1  = 150  # m/s - Assuming Constant
Cw1  = 0    # m/s
N    = 200  # rev/s
dTos_ie = 30 # K - stag temp rise for first and last stages

# Fixed Inputs
BPR = 6
mdot_a = 280
mdot_c = 40
r_c = 15
r_fan = 1.4
nf_p = 0.92
To_turb_in = 2275
nt_p = 0.92 

# Compressor Stuff
dTo_comp = 755.329 - 319.723 # K
Toi = 319.723 # K
Poi = 141400.000 # Pa 
Toe =   755.329  # K
Poe = 2121000.000  # Pa
Ca_e = Ca1

# Gas Constants
cpa = 1005 # J/kg*K
cpg = 1148 # J/kg*K
R = 287 # J/kg*K
gamma_a = 1.4
gamma_g = 4/3

# Checks
# V2/V1 >= 0.65
de_Haller_min = 0.65
# Mt <= 1.2
# phi >= 0.78
# psi <= 3.3 


# 1 - Calculate statioc at inlet
T1 = Toi - (Ca1**2) / (2*cpa)
P1 = Poi*(T1/Toi)**(gamma_a/(gamma_a-1))
rho1 = P1/(R*T1)

# Find first radii
V1_t = Mt*np.sqrt(gamma_a*R*T1)
U1_t = np.sqrt(V1_t**2 - Ca1**2)
r1_t = U1_t / (2*np.pi*N)
r1_r = r1_t *np.sqrt(1 - mdot_c/(rho1*Ca1*np.pi*r1_t**2))
r_m = (r1_t + r1_r)/2  # Constant

# Compressor Exit conditions (e represents 'exit')
Te = Toe - (Ca_e**2) / (2*cpa)
Pe = Poe*(Te/Toe)**(gamma_a/(gamma_a-1))
rho_e = Pe/(R*Te)
A_e = mdot_c / (rho_e*Ca_e)
# Assuming mean line is constant 
h_e = A_e/(2*np.pi*r_m)
re_t = r_m + h_e/2
re_r = r_m - h_e/2


# Stage number estimation
U_m = 2*np.pi*r_m*N
beta_1m = np.arctan(U_m/Ca1) # Radians
# print(np.degrees(beta_1m)) # Matches
V_1m = Ca1/np.cos(beta_1m)
V_2m = de_Haller_min*V_1m 
beta_2m = np.arccos(Ca1/(V_2m)) # Radians
# print(np.degrees(beta_2m)) # Matches

# IGNORING WorkDoneFactor
dTo_s_estimate = U_m*Ca1*(np.tan(beta_1m)-np.tan(beta_2m)) / cpa # K
stage_approx = dTo_comp / dTo_s_estimate
num_stages_approx = int(np.ceil(stage_approx))
# print('{:6.3f} => {} stages'.format(stage_approx, num_stages_approx)) Check
num_stages_approx += 1 # add one stage to account for WDF losses

# Assuming a dTos = 30K for the first and last stages
dTo_s_est_stage = dTo_comp / num_stages_approx # All stages if even dTo distribution
dTo_s_est_remaining = (dTo_comp - 2*dTos_ie) / (num_stages_approx-2) # For all intermediate stages

# =============================================================================
# # Stage by stage design
# =============================================================================
stage_num = 1
dTo = dTos_ie
nc = nc_p # Assume nc for stage = nc_poly for whole compressor:
To1 = Toi
Po1 = Poi
# def stage_calc(stage_num, U, Ca, dTo):
    # Assuming at the mean line
    
wdf = 0.98#interp_wdf(stage_num)
dCw = cpa*dTo / (wdf*U_m) # Cw2 - Cw1 = dCw
Cw2 = dCw + Cw1
beta_1 = np.arctan(U_m/Ca1)
beta_2 = np.arctan((U_m - Cw2) / Ca1)
alpha_2 = np.arctan(Cw2/Ca1)
# print(np.degrees(beta_1), np.degrees(beta_2), np.degrees(alpha_2))# Matches
# de Haller check
deHaller = np.cos(beta_1)/np.cos(beta_2)
if deHaller < de_Haller_min:
    print('de Haller check failed, V2/V1 = {:.4f}'.format(deHaller))
else: 
    print('de Haller check passed, V2/V1 = {:.4f}'.format(deHaller))
# Calculate pressuure
Po3 = Po1*(1 + nc*dTo/To1)**(gamma_a/(gamma_a - 1))
To3 = To1 + dTo
# print(To3, Po3*1e-5) # Matches
Lamda_1 = 1 - (Cw1 + Cw2)/(2*U_m) # Approximate degree of reaction


stage_num = 2
dTo = dTo_s_est_remaining
wdf = 0.93#interp_wdf(stage_num)
Lamda_2 = 0.65 # choosing based off Lamda_1 = 82.4 and Lamda_3 = 50 (why 50?)
# NOTE: Need a way to determind lamda for intermediate stages?
eq1=  (cpa*dTo) / (wdf*U_m*Ca1) # = tan(B4) - tan(B5)
eq2= 2*Lamda_2*U_m/Ca1 # = tan(B4) + tan(B5)
beta_4 = np.arctan((eq1+eq2)/2)
beta_5 = np.arctan(eq2 - np.tan(beta_4))
# print(np.degrees(beta_4)) # Matches
# print(np.degrees(beta_5)) # Matches
alpha_4 = np.arctan(U_m/Ca1  - np.tan(beta_4))
alpha_5 = np.arctan(U_m/Ca1  - np.tan(beta_5))
# print(np.degrees(alpha_4)) # Matches
# print(np.degrees(alpha_5)) # Matches
Cw4 = Ca1*np.tan(alpha_4)
Cw5 = Ca1*np.tan(alpha_5)
C4 = np.sqrt(Cw4**2 + Ca1**2)
C5 = np.sqrt(Cw5**2 + Ca1**2)

To4 = To3 # area 3 and 4 are the same
Po4 = Po3
T4 = To4 - (C4**2)/(2*cpa)
P4 = Po4*(T4/To4)**(gamma_a/(gamma_a-1))
rho4 = P4/(R*T4)

h_4 = mdot_c/(2*np.pi*rho4*Ca1*r_m)
r4_t = r_m + h_4/2
r4_r = r_m - h_4/2 

# Checking tip mach
Cw4_t = Cw4*r_m/r4_t 
C4_t = np.sqrt(Cw4_t**2 + Ca1**2)
T4_t = To4 - (C4_t**2)/(2*cpa)
V4_t = np.sqrt((2*np.pi*N*r4_t - Cw4_t)**2 + Ca1**2)
M4_t = V4_t / np.sqrt(gamma_a*R*T4_t) # Below 1.2
R4_tm = r4_t/r_m 
Lamda_2t = 1 - (1 - Lamda_2)/R4_tm**2
deHaller = np.cos(beta_4)/np.cos(beta_5)
if deHaller < de_Haller_min:
    print('de Haller check failed, V2/V1 = {:.4f}'.format(deHaller))
else: 
    print('de Haller check passed, V2/V1 = {:.4f}'.format(deHaller))
Po6 = Po4*(1 + nc*dTo/To4)**(gamma_a/(gamma_a-1)) 
To6 = To4 + dTo
