# AER E 344 Lab 7 Data Analysis
#Kyle Younge

import os
import pandas as pd
import numpy as np
# from scipy.signal import savgol_filter
# from scipy.signal import welch
# from control import mag2db
# from matplotlib import pyplot as plt

### INITS ###
DENSITY = 1.225                     # kg/m^3
VISCOSITY = 1.8305e-5               # kg/(m*s)
STD_PRESSURE = 101325               # Pa
MOTOR_SPEED = 8                     # Hz
DELTA_Y = 0.004
IN_2_MM = 25.4                      # inches to mm
MM_2_M = 0.001                      # mm ---> m
SAMPLE_FREQ = 100                  # Hz
DATA_PATH = os.path.join(os.getcwd(), 'Data')

TOTAL_PRESSURE = -2.60641           # Pa
STATIC_PRESSURE = -81.6712          # Pa
RE_TRANSITION = 5 * 10**5

# Stores all data obtained from lab setup
data_collection = []

# Stores average pressure between probes per run. Ports A and E stored in last two indices
pressure_averages = []
# Stores pressure between probes
pressure_averages_between_probes = []
# Stores velocity between probes
velocities_per_y = [] # m/s
# Store velocity in the test section
test_section_velocity = ((4.713  * MOTOR_SPEED) - 1.7961) / 2.237
print(f"U_infinity = {test_section_velocity} m/s")

# Stores integral terms for computing c_d
integral_terms = []
# Stores c_d per run using summation of integral terms
c_d = []

# Stores y_values where pressure probes are
y_values = [0.004]
y_values_mm = [4]
for i in range(37):
    y_values.append(y_values[-1] + DELTA_Y)
    y_values_mm.append(y_values_mm[-1] + 1000 * DELTA_Y)

def getAveragePressures(run_number):
    averages = []
    # Ports 1 through 16
    for i in range(1, 17):
        # Append average pressure
        averages.append(np.mean(data_collection[run_number-1].values[:,i]))

    # Ports 17 through 32
    for i in range(34, 50):
        averages.append(np.mean(data_collection[run_number-1].values[:,i]))

    # Ports 33 through 37
    for i in range(50, 54):
        averages.append(np.mean(data_collection[run_number-1].values[:,i]))

    # Return array of averages
    return averages

if True:
    # Appends average pressure data per run
    run_num = 1
    for filename in sorted(os.listdir(DATA_PATH)):
        file = os.path.join(DATA_PATH, filename)
        temp = pd.read_csv(file, header=None, encoding="utf-8")
        data_collection.append(temp - STATIC_PRESSURE)
        pressure_averages.append(getAveragePressures(run_num))
        run_num += 1

    # Computes average pressure between probes, per run
    for run in range(len(data_collection)):
        temp = [] # To be appended to pressure_between_probes
        # Per probe
        for i in range(35): # gets probe 1 to probe 37
            # print(i)
            # print(pressure_averages[run][i])
            # print(f"{pressure_averages[run][i+1]}\n\n")
            temp.append(0.5 * (pressure_averages[run][i] + pressure_averages[run][i+1]))

        pressure_averages_between_probes.append(temp)

    # Computes velocity at y_distances per run
    for run in range(len(data_collection)):
        temp = []
        for probe in range(len(pressure_averages_between_probes[run])):
            temp.append(np.sqrt(pressure_averages_between_probes[run][probe] * 2 / DENSITY))

        velocities_per_y.append(temp)

    
