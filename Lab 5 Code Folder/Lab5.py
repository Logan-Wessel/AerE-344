# AER E 344 Lab 5 Data Analysis
# Declan Green

import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

### INITS ###
CHORD_LENGTH = 0.101                # m
DENSITY = 1.225                     # kg/m^3
VISCOSITY = 1.8305e-5               # kg/(m*s)
STD_PRESSURE = 101325               # Pa
CALIBR_COEFF = 1.1                  # dimensionless
# MOTOR_SPEED = 20                    # Hz
# VELOCITY_COEFF = 1.241              # dimensionless
DATA_PATH = os.path.join(os.getcwd(), 'Data')
AIRFOIL_PATH = os.path.join(os.getcwd(), 'Airfoil Data')

# Stores all data obtained from lab setup
data_collection = []

# Angles of attack tested per run
angles_of_attack = np.deg2rad([-4, 0, 4, 6, 8, 10, 12, 14, 16]) # radians

# Stores dimensionless airfoil x-coordinate csv data
airfoil_x = pd.read_csv(os.path.join(AIRFOIL_PATH, 'x_c.csv'), header=None).values[0,:]
# Stores dimensionless airfoil y-coordinate csv data
airfoil_y = pd.read_csv(os.path.join(AIRFOIL_PATH, 'y_c.csv'), header=None).values[0,:]

# Stores average pressures per port per run. Ports A and E stored in last two indices
pressure_averages = []
# Store dynamic pressure in test section
dynamic_pressure = []
# Stores average pressure between ports
pressure_between_ports = []
# Stores measured gauge pressures between ports
gauge_pressures = []

# Stores velocity in test section per run
velocities = []
# Stores Reynolds numbers per run
reynolds_numbers = []

# Stores delta-x between ports
delta_x_i = []
# Stores delta-y between ports
delta_y_i = []
# Stores average x-coordinate between ports
avg_x_i_2 = []
# Stores average y-coordinate between ports
avg_y_i_2 = []

# Stores normal component of force between ports
normal_component_force = []
# Stores axial component of force between ports
axial_component_force = []
# Stores moment component about leading edge between ports
moment_component = []

# Stores total normal force in airfoil plane
normal_forces = []
# Stores total axial force in airfoil plane
axial_forces = []

# Stores lift force per span in x-y plane
lift_forces_per_span = []
# Stores drag force per span in x-y plane
drag_forces_per_span = []
# Stores moment about leading edge per span
moment_le_per_span = []

# Stores lift force in x-y plane
lift_forces = []
# Stores drag force in x-y plane
drag_forces = []
# Stores moment about leading edge in airfoil plane
moments_le = []

# Stores C_l per run
C_l = []
# Stores C_d per run
C_d = []
# Stores C_m per run
C_m = []
# Store C_p per run
C_p = []

### FUNCTIONS ###
# Gets average pressures per port for a specified run
def getAveragePressures(run_number):
    averages = []
    # Ports 1 through 16
    for i in range(1,17):
        # Append average pressure
        averages.append(np.mean(data_collection[run_number-1].values[:,i]))
    # Ports 17 through 32
    for j in range(34,50):
        averages.append(np.mean(data_collection[run_number-1].values[:,j]))
    # Ports 33 through 43
    for k in range(67,78):
        averages.append(np.mean(data_collection[run_number-1].values[:,k]))
    # Ports A and E as indices 43 and 44
    for n in range(78,80):
        averages.append(np.mean(data_collection[run_number-1].values[:,n]))
    # Return array of averages
    return averages

### CALCULATIONS ###
# Appends average pressure data per run
run_num = 1
for filename in sorted(os.listdir(DATA_PATH)):
    file = os.path.join(DATA_PATH, filename)
    temp = pd.read_csv(file, header=None, encoding="us-ascii")
    data_collection.append(temp)
    pressure_averages.append(getAveragePressures(run_num))
    run_num += 1

# Computes dynamic pressure in test section per run
for run in range(len(data_collection)):
    dynamic_pressure.append(CALIBR_COEFF * (pressure_averages[run][43] - pressure_averages[run][44]))

# Computes airfoil geometry data
for point in range(len(airfoil_x)):
    if point == 42:
        delta_x_i.append(airfoil_x[0] - airfoil_x[point])
        delta_y_i.append(airfoil_y[0] - airfoil_y[point])
        avg_x_i_2.append(0.5 * (airfoil_x[point] + airfoil_x[0]))
        avg_y_i_2.append(0.5 * (airfoil_y[point] + airfoil_y[0]))
    else:
        delta_x_i.append(airfoil_x[point+1] - airfoil_x[point])
        delta_y_i.append(airfoil_y[point+1] - airfoil_y[point])
        avg_x_i_2.append(0.5 * (airfoil_x[point] + airfoil_x[point+1]))
        avg_y_i_2.append(0.5 * (airfoil_y[point] + airfoil_y[point+1]))

# Computes pressure between ports per run
for run in range(len(data_collection)):
    temp = [] # To be appended to pressure_between_ports
    temp_g = [] # To be appended to gauge_pressures
    # Per airfoil port
    for i in range(len(airfoil_x)):
        if i == 42:
            temp.append(0.5 * (pressure_averages[run][i] + pressure_averages[run][0]) + STD_PRESSURE)
            temp_g.append(0.5 * (pressure_averages[run][i] + pressure_averages[run][0]))
        else:
            temp.append(0.5 * (pressure_averages[run][i] + pressure_averages[run][i+1]) + STD_PRESSURE)
            temp_g.append(0.5 * (pressure_averages[run][i] + pressure_averages[run][i+1]))

    pressure_between_ports.append(temp)
    gauge_pressures.append(temp_g)

# Computes normal, axial force component and moment component about leading edge from i-th panel
for run in range(len(data_collection)):
    temp_normal = []
    temp_axial = []
    temp_moment = []
    for point in range(len(airfoil_x)):
        temp_normal.append(pressure_between_ports[run][point] * delta_x_i[point])
        temp_axial.append(-1 * pressure_between_ports[run][point] * delta_y_i[point])
        temp_moment.append(-1 * (pressure_between_ports[run][point] * delta_x_i[point] * avg_x_i_2[point]) - (pressure_between_ports[run][point] * delta_y_i[point] * avg_y_i_2[point]))

    normal_component_force.append(temp_normal)
    axial_component_force.append(temp_axial)
    moment_component.append(temp_moment)

# Computes normal, axial force and moment about leading edge
for run in range(len(data_collection)):
    normal_forces.append(sum(normal_component_force[run]))
    axial_forces.append(sum(axial_component_force[run]))
    moment_le_per_span.append(sum(moment_component[run]))

# Computes lift and drag force per unit span, per run
for run in range(len(data_collection)):
    lift_forces_per_span.append((normal_forces[run] * np.cos(angles_of_attack[run])) - (axial_forces[run] * np.sin(angles_of_attack[run])))
    drag_forces_per_span.append((normal_forces[run] * np.sin(angles_of_attack[run])) + (axial_forces[run] * np.cos(angles_of_attack[run])))

# Computes lift and drag forces, per run
for run in range(len(data_collection)):
    lift_forces.append(lift_forces_per_span[run] * CHORD_LENGTH)
    drag_forces.append(drag_forces_per_span[run] * CHORD_LENGTH)
    moments_le.append(moment_le_per_span[run] * CHORD_LENGTH**2)

# Computes C_l, C_d, C_m, and C_p per run
for run in range(len(data_collection)):
    C_l.append(lift_forces[run] / (dynamic_pressure[run] * CHORD_LENGTH))
    C_d.append(drag_forces[run] / (dynamic_pressure[run] * CHORD_LENGTH))
    C_m.append(moments_le[run] / (dynamic_pressure[run] * CHORD_LENGTH**2))
    
# Computes C_p per run
for run in range(len(data_collection)):
    temp_pressure = []
    for point in range(len(airfoil_x)):
        temp_pressure.append((gauge_pressures[run][point] - pressure_averages[run][44]) / dynamic_pressure[run])
        
    C_p.append(temp_pressure)
    
# Computes velocities and reynolds numbers per run
for run in range(len(data_collection)):
    velocities.append(np.sqrt(dynamic_pressure[run] * 2 / DENSITY))
    reynolds_numbers.append(DENSITY * velocities[run] * CHORD_LENGTH / VISCOSITY)
    
print(velocities)
print(reynolds_numbers)

### PLOTTING ### 
# Set to True to generate plots
PLOT_FLAG = True
if PLOT_FLAG == True:
    plt.figure(1)
    plt.plot(airfoil_x, airfoil_y)
    plt.title("GA(W)-1 Airfoil Geometry")
    plt.xlabel("x/c")
    plt.ylabel("y/c")
    plt.xlim([0,1])
    plt.ylim([-.5,.5])
    plt.grid()
    plt.show()

    plt.figure(2)
    plt.plot(np.rad2deg(angles_of_attack), C_l)
    plt.title("Coefficient of Lift vs AOA")
    plt.xlabel("Angle of Attack (degrees)")
    plt.ylabel("Coefficient of Lift")
    plt.grid()
    plt.show()

    plt.figure(3)
    plt.plot(np.rad2deg(angles_of_attack), C_d)
    plt.title("Coefficient of Drag vs AOA")
    plt.xlabel("Angle of Attack (degrees)")
    plt.ylabel("Coefficient of Drag")
    plt.grid()
    plt.show()

    plt.figure(4)
    plt.plot(np.rad2deg(angles_of_attack), C_m)
    plt.title("Coefficient of Moment vs AOA")
    plt.xlabel("Angle of Attack (degrees)")
    plt.ylabel("Coefficient of Moment")
    plt.grid()
    plt.show()

    fig = plt.figure(5)
    ax = fig.add_subplot()
    ylims = [min(C_l)-0.5,max(C_l)+0.2]
    fig.subplots_adjust(right=0.75)
    twin1 = ax.twinx()
    twin2 = ax.twinx()
    # Offset the right spine of twin2.  The ticks and label have already been
    # placed on the right by twinx above.
    twin2.spines.right.set_position(("axes", 1.2))
    p1, = ax.plot(np.rad2deg(angles_of_attack), C_l, "C0", label="c_l")
    p2, = twin1.plot(np.rad2deg(angles_of_attack), C_d, "C1", label="c_d")
    p3, = twin2.plot(np.rad2deg(angles_of_attack), C_m, "C2", label="c_m")
    ax.set(ylim = ylims, xlabel="AOA (degrees)", ylabel="c_l")
    twin1.set(ylim=ylims, ylabel="c_d")
    twin2.set(ylim=ylims, ylabel="c_m")
    ax.yaxis.label.set_color(p1.get_color())
    twin1.yaxis.label.set_color(p2.get_color())
    twin2.yaxis.label.set_color(p3.get_color())
    ax.tick_params(axis='y', colors=p1.get_color())
    twin1.tick_params(axis='y', colors=p2.get_color())
    twin2.tick_params(axis='y', colors=p3.get_color())
    ax.legend(handles=[p1, p2, p3])
    plt.title("Force and Moment Coefficients vs AOA")
    plt.grid()
    plt.show()

    figure_index = 0
    for angle in angles_of_attack:
        plt.figure(figure_index + 5)
        plt.plot(airfoil_x[21:], np.negative(C_p[figure_index][21:]), color='red', label='Top Surface')
        plt.plot(airfoil_x[0:20], np.negative(C_p[figure_index][0:20]), color='blue', label='Bottom surface')
        plt.xlabel('x/c')
        plt.ylabel('-c_p')
        plt.legend()
        plt.title(f'Distribution of c_p along Airfoil at {np.rad2deg(angles_of_attack[figure_index]):.0f} degrees AOA')
        plt.grid()
        plt.show()
        figure_index += 1
