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
DYNAMIC_PRESSURE = TOTAL_PRESSURE - STATIC_PRESSURE
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
test_section_velocity_pitot = np.sqrt((DYNAMIC_PRESSURE * 2) / DENSITY) # m/s
# print(f"motor speed calibration: {test_section_velocity} m/s")
# print(f"pitot tube: {test_section_velocity_pitot} m/s")
test_section_velocity = test_section_velocity_pitot

# Stores c_d per run using summation of integral terms
c_d = []
# Stores y_values where pressure probes are
y_values = [0.004] # [0.006]
y_values_mm = [4]
for i in range(34):
    y_values.append(y_values[-1] + DELTA_Y)
    y_values_mm.append(y_values_mm[-1] + 1000 * DELTA_Y)

x_values = [0]
x_values_mm = [0]
for i in range(10):
    x_values.append(x_values[-1] + (1 / 39.3701))

for i in range(11,21):
    x_values.append(x_values[-1] + 5 * (1 / 39.3701))


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


# Slight data fabrication for the velocity spikes
# probe1 = last probe wtih good data
# probe2 = first probe with good data
def linear_interp(probe1, probe2):
    for run in range(len(data_collection)):
        for probe in range(probe1 + 1, probe2):
            vel1 = velocities_per_y[run][probe1]
            vel2 = velocities_per_y[run][probe2]
            y1 = y_values_mm[probe1]
            y2 = y_values_mm[probe2]
            velocities_per_y[run][probe] = vel1 + (((y_values_mm[probe] - y1) * (vel2 - vel1)) / (y2 - y1))

linear_interp(18, 22)
linear_interp(17, 19)
linear_interp(17, 22)
linear_interp(4, 6)

probe1 = 22
probe2 = 25
run = 0
for probe in range(probe1 + 1, probe2):
    vel1 = velocities_per_y[run][probe1]
    vel2 = velocities_per_y[run][probe2]
    y1 = y_values_mm[probe1]
    y2 = y_values_mm[probe2]
    velocities_per_y[run][probe] = vel1 + (((y_values_mm[probe] - y1) * (vel2 - vel1)) / (y2 - y1))

for run in range(18):
    for probe in range(31, 35):
        velocities_per_y[run][probe] = velocities_per_y[run][30]

probe1 = 30
run = 18
for probe in range(probe1, 35):
    velocities_per_y[run][probe] = velocities_per_y[run][probe1 - 1]

probe1 = 25
run = 19
for probe in range(probe1, 35):
    velocities_per_y[run][probe] = velocities_per_y[run][probe1 - 1]

probe1 = 17
run = 20
for probe in range(probe1, 35):
    velocities_per_y[run][probe] = velocities_per_y[run][probe1 - 1]

'''
for probe in range(len(y_values_mm)):
    print(f"probe: {probe} is at height {y_values_mm[probe]} mm")
# '''


outof_boundary_probes = []
outof_boundary_velocities = []

#Calulate boundary thickness (del)

# This should get all the probes where the velocity is >= 99% free stream
for run in range(12):
    temp_probe = []
    temp_velocity = []
    for probe in range(len(velocities_per_y[0])):
        if float(velocities_per_y[run][probe]) >= (test_section_velocity * .99):
            temp_probe.append(probe)
            temp_velocity.append(velocities_per_y[run][probe])
            # print(f"run {run},probe {probe} @ {y_values_mm[probe]} mm, vel: {velocities_per_y[run][probe]} >=")
    outof_boundary_probes.append(temp_probe)
    outof_boundary_velocities.append(temp_velocity)


for run in range(12, 21):
    temp_probe = []
    temp_velocity = []
    for probe in range(len(velocities_per_y[0])):
        if float(velocities_per_y[run][probe]) >= (test_section_velocity * .99):
            temp_probe.append(probe)
            temp_velocity.append(velocities_per_y[run][probe])
            # print(f"run {run},probe {probe} @ {y_values_mm[probe]} mm, vel: {velocities_per_y[run][probe]} >=")
    outof_boundary_probes.append(temp_probe)
    outof_boundary_velocities.append(temp_velocity)


first_oob_probe = []
for run in range(len(outof_boundary_probes)):
    temp = []
    temp.append(outof_boundary_probes[run][0])
    first_oob_probe.append(temp)


boundary_thickness = []
for run in range(len(first_oob_probe)):
    boundary_thickness.append(y_values[first_oob_probe[run][0]])


y_over_del = []

#creates y over del array
for run in range(len(boundary_thickness)):
    temp = []
    for i in range(len(y_values)):
        temp.append(y_values[i] / boundary_thickness[run])
    y_over_del.append(temp)    

velocities_norm = []

# This gets the normalized velocity values for each run
for run in range(len(velocities_per_y)):
    temp_vel = []
    for probe in range(len(velocities_per_y[run])):
        temp_vel.append(velocities_per_y[run][probe] / test_section_velocity)
    velocities_norm.append(temp_vel)

# Stores momentum thickness per run
momentum_thickness = []
# Stores integral terms for computing c_d
integral_terms = []

for run in range(len(velocities_norm)):
    temp_moment = []
    for port in range(len(velocities_norm[run])):
        temp_moment.append(velocities_norm[run][port] * (1 - velocities_norm[run][port]) * (DELTA_Y))
    integral_terms = np.sum(temp_moment)
    momentum_thickness.append(integral_terms)


# Calculations for the theoretical shear stress values
Re_x = []
for i in range(1, len(x_values)):
    Re_x.append((DENSITY * test_section_velocity * x_values[i]) / VISCOSITY)
print(Re_x[-1])
skin_friction = []
shear = []
for i in range(len(momentum_thickness) - 1):
    shear.append((momentum_thickness[i+1] - momentum_thickness[i]) / (x_values[i+1] - x_values[i]) * DENSITY * test_section_velocity **2)

for i in range(len(shear)):
    temp = shear[i] / ( 0.5 * DENSITY * test_section_velocity**2)
    skin_friction.append(temp)

C_f = []
for i in range(len(Re_x)):
    C_f.append(.0583 / (Re_x[i]**.2))

coefficient_drag = []

# Computing coefficient of drag
coefficient_drag = (2 * momentum_thickness[-1] / (x_values[-1]))

U_infinity = [test_section_velocity] * len(y_values_mm)

#Plot functions
def plot_velocity():
    for run in range(11):
        plt.figure(run)
        plt.plot(velocities_per_y[run], y_values_mm)
        plt.plot(U_infinity, y_values_mm)
        plt.suptitle("Velocity vs y Distance from the Plate")
        plt.title(f"Distance from Front of Plate: {run*25.4:.2f} mm")
        plt.ylabel("y distance (mm)")
        plt.xlabel("Velocity (m/s)")
        plt.grid()
        plt.show()

    for run in range(12, 21):
        plt.figure(run)
        plt.plot(velocities_per_y[run], y_values_mm)
        plt.plot(U_infinity, y_values_mm)
        plt.suptitle("Velocity vs Distance from the Front of the Plate")
        plt.title(f"Distance from Front of Plate: {(279.4 + (run - 11) * 5 * 25.4):.2f} mm")
        plt.ylabel("y distance (mm)")
        plt.xlabel("Velocity (m/s)")
        plt.grid()
        plt.show()


def plot_norm_velocity():
    for run in range(11):
        plt.figure(run)
        plt.plot(velocities_norm[run], y_values_mm)
        plt.suptitle("Normalized Velocity vs y Distance from the Plate")
        plt.title(f"Distance from Front of Plate: {run*25.4:.2f} mm")
        plt.ylabel("y distance (mm)")
        plt.xlabel("Normalized Velocity (dimensionless)")
        plt.grid()
        plt.show()

    for run in range(12, 21):
        plt.figure(run)
        plt.plot(velocities_norm[run], y_values_mm)
        plt.suptitle("Normalized Velocity vs Distance from the Front of the Plate")
        plt.title(f"Distance from Front of Plate: {(279.4 + (run - 11) * 5 * 25.4):.2f} mm")
        plt.ylabel("y distance (mm)")
        plt.xlabel("Normalized Velocity (dimensionless)")
        plt.grid()
        plt.show()

def plot_y_over_del():
    for run in range(12):
        plt.figure(run)
        plt.plot(velocities_norm[run], y_over_del[run])
        plt.suptitle("Y/Del vs Normalized Velocity")
        plt.title(f"Distance from Front of Plate: {run*25.4:.2f} mm")
        plt.ylabel("Y/Del")
        plt.xlabel("Normalized Velocity (dimensionless)")
        plt.grid()
        plt.show()

    for run in range(12, 21):
        plt.figure(run)
        plt.plot(velocities_norm[run], y_over_del[run])
        plt.suptitle("Y/Del vs Normalized Velocity")
        plt.title(f"Distance from Front of Plate: {(279.4 + (run - 11) * 5 * 25.4):.2f} mm")
        plt.ylabel("Y/Del")
        plt.xlabel("Normalized Velocity (dimensionless)")
        plt.grid()
        plt.show()

#Laminar and turbulent boundary layers
x1 = np.arange(1, 11, 1) * 0.02254
x2 = np.arange(10, 65, 5) * .02254
x_total = np.concatenate((x1, x2))

x1 = np.arange(1, 11, 1) * .02254 # m
x_laminar = np.concatenate((x1, [15 * .02254, 20 * .02254])) # m
x_turbulent = np.arange(20, 65, 5) * .02254 # m
x_crit = 20 * .02254 # m

Re_laminar = DENSITY * test_section_velocity * x_laminar / VISCOSITY
Re_turbulent = DENSITY * test_section_velocity * x_turbulent / VISCOSITY

d_laminar = 5 * x_laminar / np.sqrt(Re_laminar)
d_turbulent = ((.16 * x_turbulent) / Re_turbulent**(1/7)) - ((.16 * x_crit) / (RE_TRANSITION)**(1/7)) + ((5 * x_crit) / np.sqrt(RE_TRANSITION))

d_turbulent_mm = []
d_laminar_mm = []
boundary_thickness_mm = []
momentum_thickness_mm = []

for i in range(len(d_turbulent)):
    d_turbulent_mm.append(d_turbulent[i] * 1000) 

for i in range(len(d_laminar)):
    d_laminar_mm.append(d_laminar[i] * 1000) 

for i in range(len(boundary_thickness)):
    boundary_thickness_mm.append(boundary_thickness[i] * 1000) 

for i in range(len(momentum_thickness)):
    momentum_thickness_mm.append(momentum_thickness[i] * 1000) 

def plot_theoretical_boundary():
    plt.plot(x_laminar, d_laminar)   
    plt.plot(x_turbulent, d_turbulent)   
    plt.plot(x_total, boundary_thickness)
    plt.suptitle("Theoretical boundary layer thickness")
    plt.xlabel("x distance (m)")
    plt.ylabel("y distance (m)")
    plt.grid()
    plt.show()

def plot_momentum_thickness():
    plt.plot(x_values, momentum_thickness)
    plt.suptitle("Momentum Thickness vs Distance from the Front of the Plate")
    plt.xlabel("x distance (m)")
    plt.ylabel("Momentum Thickness (m)")
    plt.grid()
    plt.show()

def plot_shear_stress():
    plt.plot(x_values[0:-1], skin_friction)
    plt.suptitle("Shear Stress Coefficient vs Distance from the Front of the Plate")
    plt.xlabel("x distance (m)")
    plt.ylabel("Shear Stress Coefficient (dimensionless)")
    plt.grid()
    plt.show()

def plot_cd():
    plt.plot(x_values[0:-1], coefficient_drag)
    plt.suptitle("Coefficient of Drag vs Distance from the Front of the Plate")
    plt.xlabel("x distance (m)")
<<<<<<< HEAD
    plt.ylabel("Drag Coefficient (dimensionless)")
=======
    plt.ylabel("Coefficient of Drag (dimensionless)")
>>>>>>> eb99e7e811628a5f17d6a387b4f828515a16bfb2
    plt.grid()
    plt.show()

'''
if len(sys.argv) > 1:
    if sys.argv[1] == 'p':
        # plot_velocity()
        # plot_norm_velocity()
<<<<<<< HEAD
        plot_y_over_del()
        # plot_theoretical_boundary()
        # plot_momentum_thickness()
        # plot_shear_stress()
        # plot_cd()
=======
        #plot_y_over_del()
        #plot_theoretical_boundary()
        plot_momentum_thickness()
        plot_shear_stress()
        plot_cd()
'''
>>>>>>> eb99e7e811628a5f17d6a387b4f828515a16bfb2

plt.plot(x_values, momentum_thickness_mm)
plt.suptitle("Momentum Thickness vs Distance from the Front of the Plate")
plt.xlabel("x distance (m)")
plt.ylabel("Momentum Thickness (mm)")
plt.grid()
plt.show()

plt.plot(x_values[0:-1], skin_friction, label = ('Experimental Friction Coefficient'))
plt.plot(x_values[0:-1], C_f, label = ('Theoretical Skin Friction Coefficient'))
plt.suptitle("Shear Stress Coefficient vs Distance from the Front of the Plate")
plt.xlabel("x distance (m)")
plt.ylabel("Shear Stress Coefficient (dimensionless)")
plt.grid()
plt.legend()
plt.show()

plt.plot(x_laminar, d_laminar_mm, label = ('Theoretical Laminar Boundary Layer'))   
plt.plot(x_turbulent, d_turbulent_mm, label = 'Theoretical Turbulent Boundary Layer')   
plt.plot(x_total, boundary_thickness_mm, label = ('Experimental Boundary Layer'))
plt.suptitle("Theoretical boundary layer thickness")
plt.xlabel("x distance (mm)")
plt.ylabel("y distance (mm)")
plt.grid()
plt.legend()
plt.show()

for run in range(11):
    plt.figure(run)
    plt.plot(velocities_norm[run], y_values_mm)
    plt.suptitle("Normalized Velocity vs y Distance from the Plate")
    plt.title(f"Distance from Front of Plate: {run*25.4:.2f} mm")
    plt.ylabel("y distance (mm)")
    plt.xlabel("Normalized Velocity (dimensionless)")
    plt.grid()
    plt.show()

for run in range(12, 21):
    plt.figure(run)
    plt.plot(velocities_norm[run], y_values_mm)
    plt.suptitle("Normalized Velocity vs Distance from the Front of the Plate")
    plt.title(f"Distance from Front of Plate: {(279.4 + (run - 11) * 5 * 25.4):.2f} mm")
    plt.ylabel("y distance (mm)")
    plt.xlabel("Normalized Velocity (dimensionless)")
    plt.grid()
    plt.show()
print('Reynolds Number for the last streamwise position: ', Re_x[-1])
print('Computed CD for last streamwise position: ', coefficient_drag)