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
IN_2_MM = 25.4                      # inches to mm
MM_2_M = 0.001                      # mm ---> m
SAMPLE_FREQ = 1                  # Hz
DATA_PATH = os.path.join(os.getcwd(), 'Data')

# Throat-domain distances
throat_axis_inches = [-4.00, -1.50, -0.30, -0.18, 0.00, 0.15, 0.30, 0.45, 0.60, 0.75, 0.90, 1.05, 1.20, 1.35, 1.45] # inches
throat_axis = IN_2_MM * np.ones(len(throat_axis_inches)) * throat_axis_inches
 
# Areas per throat-domain distance
areas_inches_sq = [0.800, 0.529, 0.480, 0.478, 0.476, 0.497, 0.518, 0.539, 0.560, 0.581, 0.599, 0.616, 0.627, 0.632, 0.634] # in^2
areas = IN_2_MM**2 * np.ones(len(areas_inches_sq)) * areas_inches_sq

data = pd.read_csv('Data/Section 2.csv', header=None, encoding="utf-8")

Used_Ports = list(range(2,15, 1))
Pressures = []
for i in Used_Ports:
    Pressures.append(data.values[:,i])

velocities_points = []
# Computes velocity at y_distances per run
for time in range(len(Pressures)):
    temp = []
    velocities_points.append(np.sqrt(Pressures[time] * 2 / DENSITY))

print(velocities_points)