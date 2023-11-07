#Kyle Younge
#Lab 9 Measured

import os
import sys
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

### INITS ###
DENSITY = 1.225                     # kg/m^3
VISCOSITY = 1.8305e-5               # kg/(m*s)
STD_PRESSURE = 101325               # Pa
STD_PRES_IM = 14.6959488            # Psi
IN_2_MM = 25.4                      # inches to mm
MM_2_M = 0.001                      # mm ---> m
SAMPLE_FREQ = 1                  # Hz

GAMMA = 1.4

# Throat-domain distances
throat_axis_inches = [-4.00, -1.50, -0.30, -0.18, 0.00, 0.15, 0.30, 0.45, 0.60, 0.75, 0.90, 1.05, 1.20, 1.35, 1.45] # inches
throat_axis = IN_2_MM * np.ones(len(throat_axis_inches)) * throat_axis_inches
 
# Areas per throat-domain distance
areas_inches_sq = [0.800, 0.529, 0.480, 0.478, 0.476, 0.497, 0.518, 0.539, 0.560, 0.581, 0.599, 0.616, 0.627, 0.632, 0.634] # in^2
areas = IN_2_MM**2 * np.ones(len(areas_inches_sq)) * areas_inches_sq

data = pd.read_csv('Data/Section 2.csv', encoding="utf-8")

Used_Ports = list(range(2,17, 1))
Pressures = []
for i in Used_Ports:
    Pressures.append(data.values[:,i])

velocities_points = []
# Computes velocity at y_distances per run
for time in range(len(Pressures)):
    velocities_points.append(np.sqrt((Pressures[time] + STD_PRES_IM) * 2 / DENSITY))


def Pressure_Total_After_Throat(p_throat):
    return  p_throat * ((1.4 + 1)/2)**(1.4/.4)

def Mach_Up(pt1_by_pm):
    return np.sqrt(((pt1_by_pm**(1/3.5) - 1) / .2))

def M2_Across_Shock(M_s1):
    return (1 + (0.2 * M_s1**2)) / (1.4 * M_s1**2 - 0.2)

def Pressure_Total_Downstream(p_s2, M_s2):
    return p_s2 * ((1 + (0.2 * M_s2))**(1.4 / 0.4))

def Mach_Down(pt2_by_pm):
    return (2/0.4) * (((pt2_by_pm)**(0.4/1.4)) - 1)


print(Pressures[3,3])

'''throat = []
for i in range(len(Pressures)):
    throat.append(Pressures[3,i])'''



'''print(Pressures[0])
throat = Pressures[0]
print(throat[:,1])'''
