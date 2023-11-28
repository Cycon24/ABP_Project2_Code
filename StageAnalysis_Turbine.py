"""
@author: ZainHita
"""

# Imports ----------
import numpy as np


# --------------------------------------------------
# --------------------------------------------------


# Gas Constants ----------
cpa     = 1.005 * 1e3 # J/kg*K
cpg     = 1.148 * 1e3 # J/kg*K
R       = 287 # J/kg*K
gamma_a = 1.4
gamma_g = 4/3


# --------------------------------------------------
# --------------------------------------------------


# Fixed Inputs ----------

# From Fig 1.
BPR    = 6
mdot_a = 280 # kg/s

# From Fig 2.
r_c   = 15.0
r_fan = 1.4
nf_p  = 0.92
nt_p  = 0.92

# From Pg 1.
Po13 = Po2 = Po1 = Pa = 1.01 # bar 
To13 = To2 = To1 = Ta = 288 # K

# Chosen/Assumed
nc_p = 0.90
C1 = Ca1 = 150 # m/s, no inlet swirl assumed.
N  = 200 # rev/s
C3 = C2 = 150 # m/s
Lambda = 0.98
dTos_ie = 30 # K

# Other Fixed Inputs
dTos_comp = 435.606 # K
Toi = 319.723 # K
Poi = 141400.000 # Pa 
Toe =   755.329  # K
Poe = 2121000.000  # Pa
Ca_e = Ca1
To_turb_in = 2275 # K

# --------------------------------------------------
# --------------------------------------------------


# Variating Inputs ----------

# From Fig 2.
V2_over_V1 = 0.65 # >= 0.65 , de Haller Criteria
nc_p       = 0.90 # >= 0.90
Mt         = 1.15 # <= 1.20


# --------------------------------------------------
# --------------------------------------------------


# Sample Calc ----------

# From Cycle Level Calcs
Po13           = r_fan * Po1 # bar
To13_minus_To2 = To2 * (r_fan**((gamma_a - 1)/(gamma_a*nf_p)) - 1) # K
To13           = To2 + To13_minus_To2 # K
mdot_c         = mdot_a / (BPR + 1) # kg/s
Po3            = Po13 * r_fan # bar
To3            = To13 * r_fan**((gamma_a - 1)/(gamma_a*nc_p)) # K

# From Compressor Calcs
T1   = To1 - Ca1**2 / (2*cpa) # K
P1   = Po1 * (T1/To1)**(gamma_a/(gamma_a - 1)) # bar
rho1 = P1 * 1e5 / (R*T1) # kg/m^3
V1t  = Mt * np.sqrt(gamma_a*R*T1) # m/s
Ut   = np.sqrt(V1t**2 + Ca1**2) # m/s
rt   = Ut / (2*np.pi*N) # m
rr   = rt * np.sqrt(1 - mdot_c / (rt**2 * rho1 * Ca1 * np.pi)) # m
rm   = (rt + rr) / 2 # m
T3   = 744.13499 # K
P3   = 20.13006 # bar
rho3 = P3 / (R*T3) # kg/m^3
A3   = 0.02829 # m^2
h    = 0.019135 # m
rte  = rm + h/2 # m
rre  = rm - h/2 # m
# Estimating number of Stages
Um   = 2*np.pi*rm*N
# Relative Velocity
Ca     = Ca1
beta1  = np.arctan(Um/Ca) # rad
V1     = Ca/np.cos(beta1) # m/s
V2     = V2_over_V1*V1 # m/s
beta2  = np.arccos(Ca/V2) # rad
dTos   = Um*Ca*(np.tan(beta1) - np.tan(beta2)) / cpa # K, work factor ingored
stages = int(np.ceil(dTos_comp / dTos)) + 1 # 1 Stage added to account for wdf # Came out to be 11, should have been 12?
dTos_per_stage = dTos_comp / stages # K
dTos_per_rem   = (dTos_comp - 2*dTos_ie) / (stages - 2) # K
# Stage by Stage Design
stage_num = 1
U_m = Um
Cw1 = 0
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
if deHaller < V2_over_V1:
    print('de Haller check failed, V2/V1 = {:.4f}'.format(deHaller))
else: 
    print('de Haller check passed, V2/V1 = {:.4f}'.format(deHaller))
# Calculate pressuure
Po3 = Po1*(1 + nc*dTo/To1)**(gamma_a/(gamma_a - 1))
To3 = To1 + dTo
# print(To3, Po3*1e-5) # Matches
Lamda_1 = 1 - (Cw1 + Cw2)/(2*U_m) # Approximate degree of reaction


stage_num = 2
dTo = dTos_per_rem
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

r_m = rm

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
if deHaller < V2_over_V1:
    print('de Haller check failed, V2/V1 = {:.4f}'.format(deHaller))
else: 
    print('de Haller check passed, V2/V1 = {:.4f}'.format(deHaller))
Po6 = Po4*(1 + nc*dTo/To4)**(gamma_a/(gamma_a-1)) 
To6 = To4 + dTo

###################################################################
"""FOR TURBINE"""
###################################################################

# Gas Constants ----------
cpa     = 1.005 * 1e3 # J/kg*K
cpg     = 1.148 * 1e3 # J/kg*K
R       = 287 # J/kg*K
gamma_a = 1.4
gamma_g = 4/3

# Required Values
dTos_per_stage = 182.33 # K
To_turb_in = 2275 # K=

# Thus, presure Ratio for per stage
Po4cycle = 21.21 # bar
Po3_over_Po4cycle = (1 - dTos_per_stage/To_turb_in)**(gamma_g/(nt_p*(gamma_g - 1)))
Po3 = Po4cycle * Po3_over_Po4cycle

# Other Assumed values
Ca1     = 150 # m/s
psi     = 3.3
phi     = 0.78
deg_rxn = 0.5

# Velocities at Mean
Um = np.sqrt(2*cpg*dTos_per_stage/psi)
Ca = phi * Um

# Gas angle and Swirl Angles Respectively
beta2  = np.arctan(0.5/phi * (psi/2 - 2*deg_rxn)) # rad
beta3  = np.arctan(0.5/phi * (psi/2 + 2*deg_rxn)) # rad
alpha2 = np.arctan(np.tan(beta2) + 1/phi) # rad
alpha3 = np.arctan(np.tan(beta3) - 1/phi) # rad

# Mean Radius
rm = Um / (2*np.pi*N) # m

# Obtaining Bade Height
C3          = Ca/np.cos(alpha3) # m/s
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

## THE END.

