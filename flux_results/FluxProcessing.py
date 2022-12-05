# -*- coding: utf-8 -*-
"""
Created on Wed May 18 21:15:04 2022

@author: Felipe
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import gridspec
from scipy.interpolate import griddata

'''
This module contains functions to do flux calculation using the Tributary Area Method
Please note the flux_integral function is customized for a plane with normal (1,0,0)
'''

# Function to calculate equivalent area of each grid element
def dA(grid_y,grid_z):
    min_y = np.min(grid_y[:,0])
    max_y = np.max(grid_y[:,0])
    min_z = np.min(grid_z.T[:,0])
    max_z = np.max(grid_z.T[:,0])
    length_y = max_y - min_y    
    length_z = max_z - min_z  
    dAreas = np.zeros((np.shape(grid_y)[0],np.shape(grid_z)[1]))
    dA_unit = length_y/(np.size(dAreas[0,:]) - 1) * length_z/(np.size(dAreas[:,0]) - 1)
    for i in range(0,np.shape(grid_y)[0]): # Columns (y - axis)
        for j in range(0,np.shape(grid_z)[0]): # Rows (z - axis)
            y = grid_y[i,0]
            z = grid_z.T[j,0]
            if ((y == min_y) or (y == max_y)) and ((z == min_z) or (z == max_z)):
                dAreas[i,j] = dA_unit/4
            elif (y == min_y) or (y == max_y):
                dAreas[i,j] = dA_unit/2
            elif (z == min_z) or (z == max_z):
                dAreas[i,j] = dA_unit/2
            else:
                dAreas[i,j] = dA_unit       
    return dAreas

# This function calculates the area integral for a plane with normal (1, 0, 0)
# This could be updated for a generic plane
def flux_integral(grid_y,grid_z,flux,dAreas): 
    if type(dAreas) == int:
        dAreas = dA(grid_y,grid_z)
    flux_int = 0
    for i in range(0,np.shape(grid_y)[0]): # Columns (y - axis)
        for j in range(0,np.shape(grid_y)[0]): # Rows (z - axis)
            if np.isnan(flux[i,j]):
                flux_int += 0
            else:
                flux_int += flux[i,j]*dAreas[i,j]
    return flux_int


    