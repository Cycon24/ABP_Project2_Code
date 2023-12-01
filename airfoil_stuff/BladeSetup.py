# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 11:32:40 2023

@author: cycon
"""
from AirfoilTransformer import AirfoilTransformer
import numpy as np
import matplotlib.pyplot as plt
# Compressor Outputs
r_r = 0.210 # m
theta_1r = 43.819 # deg
theta_2r = -5.897 # deg

r_m = 0.29 # m
theta_1m = 45.292 # deg
theta_2m = 0 # deg

r_t = 0.370 # m
theta_1t = 46.679 # deg
theta_2t = 6.397  # deg


# Params
c = 0.0533 # m - chord length
s = 0.9*c  # m - blade pitch (s/c = 0.9)
n = 2*np.pi*r_m/s
print(f'Number of blades: {n} = {np.ceil(n)}')

airfoil_file_location = "C:\\Users\\cycon\\Documents\\ABP_P2\\airfoil_stuff\\airfoils\\"
airfoil_file_name = "NACA-9509"
camberline_file_loc = "C:\\Users\\cycon\\Documents\\ABP_P2\\airfoil_stuff\\airfoils\\Camberline.txt"

# Initiate each airfoil
root = AirfoilTransformer(airfoil_file_location, airfoil_file_name)
mean = AirfoilTransformer(airfoil_file_location, airfoil_file_name)
tip  = AirfoilTransformer(airfoil_file_location, airfoil_file_name)

LE_camber_deg, TE_camber_deg = root.camberline_angles(camberline_file_loc)

root_rot = -(theta_1r + theta_2r) #- LE_camber_deg
mean_rot = -(theta_1m + theta_2m) #- LE_camber_deg
tip_rot  = -(theta_1t + theta_2t) # - LE_camber_deg

scaleFac = 1
x_transFac = -scaleFac/2
r_r = r_r / c
r_m = r_m / c
r_t = r_t / c

root.Transform(scale=scaleFac, translations=[x_transFac,0,r_r], rotations=[0,0,root_rot], rotateFirst=False)
mean.Transform(scale=scaleFac, translations=[x_transFac,0,r_m], rotations=[0,0,mean_rot], rotateFirst=False)
tip.Transform(scale=scaleFac, translations=[x_transFac,0,r_t], rotations=[0,0,tip_rot], rotateFirst=False)

fig = plt.figure("Airfoil Transformer")
ax = fig.add_subplot(121, projection='3d')

root.PlotAirfoils(fig=fig, ax=ax)
tip.PlotAirfoils(fig=fig, ax=ax)
mean.PlotAirfoils(fig=fig, ax=ax)

root.SaveAsTXT("root")
tip.SaveAsTXT("tip")
mean.SaveAsTXT("mean")
