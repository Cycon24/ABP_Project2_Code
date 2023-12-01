# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 11:32:40 2023

@author: cycon
"""
from AirfoilTransformer import AirfoilTransformer
import numpy as np
import matplotlib.pyplot as plt
from thetaCalc import CalculateThetas



# Compressor Outputs
r_r = 0.222 # m
beta_1r = 42.246
beta_2r = 10.628

r_m = 0.265 # m
beta_1m = 49.912
beta_2m = 30.262

r_t = 0.307 # m
beta_1t = 55.491
beta_2t = 43.039

# Params
AR = 3
s_c = 0.9 # s/c - pitch to chord ratio
h = r_t - r_r 
c = h/AR   # m - chord length
s = s_c*c  # m - blade pitch (s/c = 0.9)

theta_1r, theta_2r = CalculateThetas(beta_1r, beta_2r, h, s_c)
theta_1m, theta_2m = CalculateThetas(beta_1m, beta_2m, h, 0.8)
theta_1t, theta_2t = CalculateThetas(beta_1t, beta_2t, h, 0.6)


# Desired Turn Angle - Phi (deg)
phi_r = theta_1r - theta_2r
phi_m = theta_1m - theta_2m
phi_t = theta_1t - theta_2t
print("Phi (φ)")
print(" root: {:.3f}".format(phi_r))
print(" mean: {:.3f}".format(phi_m))
print(" tip:  {:.3f}".format(phi_t))


n = 2*np.pi*r_m/s
print(f'Number of blades: {n} = {np.ceil(n)}')

airfoil_file_location = "C:\\Users\\cycon\\Documents\\ABP_P2\\airfoil_stuff\\airfoils\\"
airfoil_name_r = "NACA-8509"
airfoil_name_m = "NACA-0009"
airfoil_name_t = "NACA-2509"
camberline_file_tag = "_Camberline.txt"
camber_r_file = airfoil_file_location + airfoil_name_r + camberline_file_tag
camber_m_file = airfoil_file_location + airfoil_name_m + camberline_file_tag
camber_t_file = airfoil_file_location + airfoil_name_t + camberline_file_tag

# Initiate each airfoil
root = AirfoilTransformer(airfoil_file_location, airfoil_name_r)
mean = AirfoilTransformer(airfoil_file_location, airfoil_name_m)
tip  = AirfoilTransformer(airfoil_file_location, airfoil_name_t)

# Get Airfoil camber angles - Alpha (deg)
alpha_LE_r, alpha_TE_r = root.camberline_angles(camber_r_file, 10)
alpha_LE_m, alpha_TE_m = mean.camberline_angles(camber_m_file, 10)
alpha_LE_t, alpha_TE_t = tip.camberline_angles(camber_t_file, 10)

# Find Airfoil Turn Angle - Lamda
lamda_r = alpha_LE_r - alpha_TE_r
lamda_m = alpha_LE_m - alpha_TE_m
lamda_t = alpha_LE_t - alpha_TE_t
print("Root - {}:\n a_LE: {:.2f}\n a_TE: {:.2f}\n λ:    {:.2f}".format(airfoil_name_r, alpha_LE_r, alpha_TE_r, lamda_r))
print("Mean - {}:\n a_LE: {:.2f}\n a_TE: {:.2f}\n λ:    {:.2f}".format(airfoil_name_m, alpha_LE_m, alpha_TE_m, lamda_m))
print("Tip - {}:\n a_LE: {:.2f}\n a_TE: {:.2f}\n λ:    {:.2f}".format(airfoil_name_t, alpha_LE_t, alpha_TE_t, lamda_t))

# Check error between Phi and Lamda
phi_lam_error_r = phi_r - lamda_r
phi_lam_error_m = phi_m - lamda_m
phi_lam_error_t = phi_t - lamda_t
print("Errors (φ-λ):")
print(" e_r = {:.6f}".format(phi_lam_error_r))
print(" e_m = {:.6f}".format(phi_lam_error_m))
print(" e_t = {:.6f}".format(phi_lam_error_t))

# Chord Line angle - Psi
psi_r = theta_1r - alpha_LE_r
psi_m = theta_1m - alpha_LE_m
psi_t = theta_1t - alpha_LE_t
print("Psi (ψ)")
print(" root: {:.3f}".format(psi_r))
print(" mean: {:.3f}".format(psi_m))
print(" tip:  {:.3f}".format(psi_t))



scaleFac = 1
x_transFac = 0 #scaleFac
r_r = r_r / c
r_m = r_m / c
r_t = r_t / c

root.Transform(scale=scaleFac, translations=[x_transFac,0,r_r], rotations=[0,0,psi_r], rotateFirst=False)
mean.Transform(scale=scaleFac, translations=[x_transFac,0,r_m], rotations=[0,0,psi_m], rotateFirst=False)
tip.Transform(scale=scaleFac, translations=[x_transFac,0,r_t], rotations=[0,0,psi_t], rotateFirst=False)

fig = plt.figure("Airfoil Transformer")
ax = fig.add_subplot(121, projection='3d')

root.PlotAirfoils(fig=fig, ax=ax, LabelTag="Root")
mean.PlotAirfoils(fig=fig, ax=ax, LabelTag="Mean")
tip.PlotAirfoils(fig=fig, ax=ax, LabelTag="Tip")

root.SaveAsTXT("root_stage8")
tip.SaveAsTXT("tip_stage8")
mean.SaveAsTXT("mean_stage8")
