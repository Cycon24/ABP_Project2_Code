# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 16:04:56 2023

@author: cycon
"""
import numpy as np



def CalculateThetas(beta1, beta2, h, s_c, maxCamberLoc=0.5, AR=3):
    sigma = beta1 - beta2
    
    c = h/AR
    
    # Point of max camber (Chosen based on airfoil or vice versa)
    a = maxCamberLoc*c
    m = 0.23*(2*a/c)**2 + 0.1*(beta2/50) # beta in degrees
    
    theta1 = sigma/(1-m*np.sqrt(s_c))
    delta = m*np.sqrt(s_c)*theta1
    
    theta2 = beta2 - delta
    
    return theta1, theta2

if __name__ == "__main__":
    # Test against nathan's hand calcs
    h = 0.271-0.258
    beta1 = 44.767 # degrees
    beta2 = 18.669  # degrees 
    s_c = 0.9


    print(CalculateThetas(beta1, beta2, h, s_c))


