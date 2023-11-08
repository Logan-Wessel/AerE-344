#Kyle Younge
#Lab 9 Measured

import os
import sys
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve

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
throat_axis_inches = [-4, -1.50, -0.30, -0.18, 0.00, 0.15, 0.30, 0.45, 0.60, 0.75, 0.90, 1.05, 1.20, 1.35, 1.45] # inches
throat_axis = IN_2_MM * np.ones(len(throat_axis_inches)) * throat_axis_inches
# Areas per throat-domain distance
areas_inches_sq = [0.800, 0.529, 0.480, 0.478, 0.476, 0.497, 0.518, 0.539, 0.560, 0.581, 0.599, 0.616, 0.627, 0.632, 0.634] # in^2
areas = IN_2_MM**2 * np.ones(len(areas_inches_sq)) * areas_inches_sq


DATA_PATH = os.path.join(os.getcwd(), 'Data')
file = os.path.join(DATA_PATH, 'Section 2.csv')
data = pd.read_csv(file, sep=',', encoding="utf-8", header = None)


# Put Pressure data from the ports we're using into an array
Used_Ports = range(3,18, 1)
Gauge_Pressures = []


for i in Used_Ports:
   Gauge_Pressures.append(data.values[:,i])


# Convert to Absolute Pressure in pascals
Pressures = []
for i in range(len(Gauge_Pressures)):
   temp = []
   for j in range(len(Gauge_Pressures[i])):
      temp.append(((Gauge_Pressures[i][j]) * 6894.76) + 101325)
   Pressures.append(temp) 

# Define functions for Total Pressure and Mach Numbers along the nozzle

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

def sp_inverse(sp, g):
    '''
    Parameters:
    sp: Pressure ratio p0/p
    g: Specific heat ratio

    Returns:
    m: Mach number
    '''
    sp_residual = lambda m, sp_term, g : sp_term - st(m, g)**(g/(g-1))
    m = fsolve(sp_residual, x0=4, args=(sp, g))[0]
    return m
st = lambda m, g : 1 + ((g - 1)/2) * m**2

# Pull Throat Pressures Column from Array
Pressures_Throat = Pressures[0]

# Compute Total Pressures Upstream from Shock
Pressures_Total_Up = []
for i in range(len(Pressures_Throat)):
   Pressures_Total_Up.append(Pressure_Total_After_Throat(Pressures_Throat[i]))


# Compute Mach Upstream for all Data Points
Mach_Upstream = []
for i in range(len(Pressures)):
   temp = []
   for j in range(len(Pressures[i])):
      temp.append(Mach_Up(Pressures_Total_Up[j])/(Pressures[i][j]))
   Mach_Upstream.append(temp)


# Compute M2 across shocks at every point
M2_Shock = []
for i in range(len(Mach_Upstream)):
   temp = []
   for j in range(len(Mach_Upstream[i])):
      temp.append(M2_Across_Shock(Mach_Upstream[i][j]))
   M2_Shock.append(temp)


for i in range(len(Pressures)):
   # print(len(Pressures[i]))
   # print(Pressures[i])
   # print()
   for j in range(len(Pressures[i])):
       print(Pressures[i][j])


frames = []
for i in range(len(Pressures[1])):
   frames.append([item[i] for item in Pressures])


# Case 1 - Frame 23
# Case1_TP


#total_pressure = the_total_pressure_you_want * np.ones(len(size of the thing you want))



# Case 2 - Frame 68
# Total Pressure Before Shock
Case2_TP = []
for i in range(len(frames[68])):
   Case2_TP.append(Pressure_Total_After_Throat(Pressures_Throat[68]))

#Mach before shock
Case2_M = []
for i in range(len(frames[68])):
   Case2_M.append(sp_inverse(Pressures_Throat[68] / frames[68][i], 1.4))


# Case 3 - Frame 128



# Case 4 - Frame 161
Case4_TP = []
for i in range(14):
   Case4_TP.append(Pressure_Total_After_Throat(Pressures_Throat[160]))

Case4_M = []
for i in range(14):
   Case4_M.append(sp_inverse(Pressures_Throat[161] / frames[161][i], 1.4))

Case4_M.append(M2_Across_Shock(Case4_M[13]))

Case4_TP.append(Pressure_Total_Downstream((Case4_TP[13]), Case4_M[14])) 

# Case 5 - Frame 227
Case5_TP = []
for i in range(len(frames[68])):
   Case5_TP.append(Pressure_Total_After_Throat(frames[68][i]))


# Case 6 - Frame 258
Case6_TP = []
Case6_TP.append(Pressure_Total_After_Throat(Pressures_Throat[258]))

Case6_M = []

Case6_M.append(sp_inverse(Pressures_Throat[258] / frames[258][i], 1.4))

Case6_M.append(M2_Across_Shock(Case4_M[12]))

for i in range (13):
   Case6_M.append(Mach_Down(frames[258][2] / frames[258][i+2]))

for i in range(14):
   Case6_TP.append(Pressure_Total_Downstream(frames[258][2], Case6_M[2]))


'''
frame_num = 227
plt.figure(1)
plt.plot(throat_axis, frames[frame_num], label='pressures per port in frame {frame_num}')
plt.legend()
plt.grid()
plt.grid()
plt.show()
plt.close()
'''

plt.figure(1)
plt.plot(throat_axis, Case2_TP, label='CV3 - Case 2')
plt.plot(throat_axis, Case4_TP, label='CV2 - Case 4')
plt.plot(throat_axis, Case6_TP, label='CV1 - Case 6')
plt.suptitle('Total Pressure vs Distance Along Nozzle')
plt.xlabel('Distance Along Nozzle (mm)')
plt.ylabel('Total Pressure (Pa)')
plt.grid()
plt.legend()
plt.show()

plt.figure(2)
plt.plot(throat_axis, Case2_M, label='CV3 - Case 2')
plt.plot(throat_axis, Case4_M, label='CV2 - Case 4')
plt.plot(throat_axis, Case6_M, label='CV1 - Case 6')
plt.suptitle('Mach vs Distance Along Nozzle')
plt.xlabel('Distance Along Nozzle (mm)')
plt.ylabel('Mach Number')
plt.grid()
plt.legend()
plt.show()

plt.figure(3)
plt.plot(throat_axis, frames[23], label='Underexpanded - Case 1')
plt.plot(throat_axis, frames[68], label='CV3 - Case 2')
plt.plot(throat_axis, frames[128], label='Overexpanded - Case 3')
plt.plot(throat_axis, frames[161], label='CV3 - Case 4')
plt.plot(throat_axis, frames[227], label='Shock in Nozzle - Case 5')
plt.plot(throat_axis, frames[258], label='CV1 - Case 6')
plt.suptitle('Static Pressure vs Distance Along Nozzle')
plt.xlabel('Distance Along Nozzle (mm)')
plt.ylabel('Total Pressure (Pa)')
plt.grid()
plt.legend()
plt.show()