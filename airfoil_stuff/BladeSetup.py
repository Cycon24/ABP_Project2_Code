# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 11:32:40 2023

@author: cycon
"""
from AirfoilTransformer import AirfoilTransformer
import numpy as np
# Compressor Outputs
r_m = 0.218

# Params
c = 0.0533 # m - chord length
s = 0.9*c  # m - blade pitch (s/c = 0.9)
n = 2*np.pi

airfoil_file_location = "C:\\Users\\cycon\\Documents\\ABP_P2\\airfoil_stuff\\airfoils\\"
airfoil_file_name = "NACA-9509"
camberline_file_loc = "C:\\Users\\cycon\\Documents\\ABP_P2\\airfoil_stuff\\airfoils\\Camberline.txt"


root = AirfoilTransformer(airfoil_file_location, airfoil_file_name)
mean = AirfoilTransformer(airfoil_file_location, airfoil_file_name)
tip  = AirfoilTransformer(airfoil_file_location, airfoil_file_name)

root.PlotAirfoils()
root.camberline_angles(camberline_file_loc)