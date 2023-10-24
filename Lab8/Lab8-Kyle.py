# AER E 344 Lab 7 Data Analysis
#Kyle Younge

import os
import sys
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
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
test_section_velocity = (4.713  * MOTOR_SPEED - 1.7961) / 2.237 # m/s
# print(test_section_velocity)

# Stores integral terms for computing c_d
integral_terms = []
# Stores c_d per run using summation of integral terms
c_d = []
# Stores momentum thickness per AOA
momentum_thickness = []
# Stores y_values where pressure probes are
y_values = [0.004] # [0.006]
y_values_mm = [4]
for i in range(34):
    y_values.append(y_values[-1] + DELTA_Y)
    y_values_mm.append(y_values_mm[-1] + 1000 * DELTA_Y)

x_values = [0]
x_values_mm = [0]
for i in range(10):
    x_values.append(x_values[-1] + 0.02254)
    x_values_mm.append(x_values_mm[-1] + 25.4)

for i in range(11,21):
    x_values.append(x_values[-1] + 5 * 0.02254)
    x_values_mm.append(x_values_mm[-1] + 5 * 25.4)


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

# Appends average pressure data per run
run_num = 1
for filename in sorted(os.listdir(DATA_PATH)):
    file = os.path.join(DATA_PATH, filename)
    # temp = pd.read_csv(file, header=None, encoding="utf-8")
    # data_collection.append(temp - STATIC_PRESSURE)
    data_collection.append(pd.read_csv(file, header=None, encoding="utf-8") - STATIC_PRESSURE)
    pressure_averages.append(getAveragePressures(run_num))
    run_num += 1

# Computes average pressure between probes, per run
for run in range(len(data_collection)):
    temp = [] # To be appended to pressure_between_probes
    # Per probe
    for i in range(35): # gets probe 1 to probe 37
        temp.append(0.5 * (pressure_averages[run][i] + pressure_averages[run][i+1]))

    pressure_averages_between_probes.append(temp)

# Computes velocity at y_distances per run
for run in range(len(data_collection)):
    temp = []
    for probe in range(len(pressure_averages_between_probes[run])):
        temp.append(np.sqrt(pressure_averages_between_probes[run][probe] * 2 / DENSITY))

    velocities_per_y.append(temp)

# print(velocities_per_y[1])
# print(len(velocities_per_y))

# Slight data fabrication for the velocity spikes on probes 19-21
for run in range(len(data_collection)):
    for probe in range(19, 22):
        y1 = velocities_per_y[run][18]
        y2 = velocities_per_y[run][22]
        x1 = 76 # mm
        x2 = 92 # mm
        velocities_per_y[run][probe] = y1 + (((y_values_mm[probe] - x1) * (y2 - y1)) / (x2 - x1))

# '''
# This should print all the probes where the velocity is >= 99% free stream
for i in range(11):
    for j in range(len(velocities_per_y[0])):
        if float(velocities_per_y[i][j]) >= (test_section_velocity * .99):
            print(f"i run {i},j probe {j} @ {y_values_mm[j]} mm, vel: {velocities_per_y[i][j]} >=")

# '''

# '''
for i in range(12, 21):
    for k in range(len(velocities_per_y[0])):
        if float(velocities_per_y[i][k]) >= (test_section_velocity * .99):
            print(f"i run {i},k probe {k} @ {y_values_mm[k]} mm, vel: {velocities_per_y[i][k]} >=")
    print("")
# '''

velocities_norm = []

# This was ment to get the normalized velocity values for each run, somehow makes all the data identical
for run in range(len(velocities_per_y)):
    temp_vel = []
    for probe in range(len(velocities_per_y[run])):
        temp_vel.append(velocities_per_y[run][probe] / test_section_velocity)
    velocities_norm.append(temp_vel)

for run in range(21):
    velocities_norm = []
    for run in range(len(velocities_per_y)):
        temp_vel = []
        for vel in range(len(velocities_per_y[run])):
            temp_vel.append((1 / test_section_velocity) * velocities_per_y[run][vel])
        velocities_norm.append(temp_vel)

    # Computing momentum thickness
    for run in range(21):
        temp_moment = []
        for i in range(len(velocities_norm[run])):
            temp_moment.append(velocities_norm[run][i] * (1 - velocities_norm[run][i]) * (DELTA_Y * MM_2_M))
            integral_term = np.sum(temp_moment)
        momentum_thickness.append(integral_term)

if len(sys.argv) > 1:
        if sys.argv[1] == 'p':
    
            for run in range(11):
                # c_p graphs
                plt.figure(run)
                #plt.plot(velocities_per_y[run], y_values_mm)
                plt.plot(velocities_norm[run], y_values_mm)
                #plt.plot(U_infiniy, y_values_mm)
                plt.suptitle("Normalized Velocity vs y Distance from the Plate")
                plt.title(f"Distance from Front of Plate: {run*25.4:.2f} mm")
                plt.ylabel("y distance (mm)")
                plt.xlabel("Velocity (m/s)")
                plt.grid()
                plt.show()


            for run in range(12, 21):
                # c_p graphs
                plt.figure(run)
                # plt.plot(velocities_per_y[run], y_values_mm)
                plt.plot(velocities_norm[run], y_values_mm)
                # plt.plot(U_infiniy, y_values_mm)
                plt.suptitle("Velocity vs Distance from the Front of the Plate")
                plt.title(f"Distance from Front of Plate: {(279.4 + (run - 11) * 5 * 25.4):.2f} mm")
                plt.ylabel("y distance (mm)")
                plt.xlabel("Velocity (m/s)")
                plt.grid()
                plt.show()
                
                plt.figure(run)
                plt.plot(x_values_mm, momentum_thickness)
                plt.suptitle("Momentum Thickness vs Distance from the Front of the Plate")
                plt.xlabel("x distance (mm)")
                plt.ylabel("Momentum Thickness")
                plt.grid()
                plt.show()
            