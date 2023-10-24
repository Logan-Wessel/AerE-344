# AER E 344 Lab 6 Data Analysis
# Declan Green

import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
# from sklearn.metrics import r2_score

### INITS ###
DENSITY = 1.225                     # kg/m^3
VISCOSITY = 1.8305e-5               # kg/(m*s)
STD_PRESSURE = 101325               # Pa
CALIBR_COEFF = 1.1                  # dimensionless
VELOCITY_COEFF = 1.241              # dimensionless
MOTOR_SPEED = 15                    # Hz
DELTA_Y =  0.002                    # meters
CHORD_LENGTH = 0.101                # meters
DATA_PATH = os.path.join(os.getcwd(), 'AER E 344/Lab 6/Data')

# Stores all data obtained from lab setup
data_collection = []

# Angles of attack tested per run
angles_of_attack = np.deg2rad([-4, 0, 4, 8, 12]) # radians

# Stores average pressure between probes per run. Ports A and E stored in last two indices
pressure_averages = []
# Stores pressure between probes
pressure_averages_between_probes = []
# Stores velocity between probes
velocities_per_y = [] # m/s
# Store velocity in the test section
test_section_velocity = VELOCITY_COEFF * MOTOR_SPEED

# Store c_p values between each probe
c_p = []

# Stores integral terms for computing c_d
integral_terms = []
# Stores c_d per run using summation of integral terms
c_d = []

# Stores y_values where pressure probes are
y_values = [-0.046]
y_values_mm = [-46]
for i in range(45):
    y_values.append(y_values[-1] + DELTA_Y)
    y_values_mm.append(y_values_mm[-1] + 1000 * DELTA_Y)

lab_5_c_d_data = [0.028494216717464306, 0.013006934230472876, -0.10890904323172532, -0.07320476145990733, 0.104224774953789]
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
    for k in range(67,81):
        averages.append(np.mean(data_collection[run_number-1].values[:,k]))
    # Ports A and E as indices 46 and 47
    for n in range(81,83):
        averages.append(np.mean(data_collection[run_number-1].values[:,n]))
    # Return array of averages
    return averages

### CALCULATIONS ###
runFlag = False
if runFlag == True:
    # Appends average pressure data per run
    run_num = 1
    for filename in sorted(os.listdir(DATA_PATH)):
        file = os.path.join(DATA_PATH, filename)
        temp = pd.read_csv(file, header=None, encoding="utf-8")
        data_collection.append(temp)
        pressure_averages.append(getAveragePressures(run_num))
        run_num += 1

    # Computes average pressure between probes, per run
    for run in range(len(data_collection)):
        temp = [] # To be appended to pressure_between_probes
        # Per probe
        for i in range(46): # gets probe 1 to probe 46
            temp.append(0.5 * (pressure_averages[run][i] + pressure_averages[run][i+1]))

        pressure_averages_between_probes.append(temp)

    # Computes c_p per space between probes
    for run in range(len(data_collection)):
        temp = [] # To be appended to c_p
        for probe in range(len(pressure_averages_between_probes[run])):
            temp.append((pressure_averages_between_probes[run][probe] - pressure_averages[run][47])/(CALIBR_COEFF * (pressure_averages[run][46] - pressure_averages[run][47])))

        c_p.append(temp)
        
    # Computes velocity at y_distances per run
    for run in range(len(data_collection)):
        temp = []
        for probe in range(len(pressure_averages_between_probes[run])):
            temp.append(np.sqrt(pressure_averages_between_probes[run][probe] * 2 / DENSITY))
            
        velocities_per_y.append(temp)

    # Computes c_d per run
    for run in range(len(data_collection)):
        temp = [] # To be appended to integral_terms
        for probe in range(len(pressure_averages_between_probes[run])):
            temp.append(((velocities_per_y[run][probe] / test_section_velocity) * (1 - velocities_per_y[run][probe] / test_section_velocity)) * DELTA_Y)

        integral_terms.append(temp)
        c_d.append((2 / CHORD_LENGTH) * sum(integral_terms[run]))
        
    for run in range(len(data_collection)):
        # c_p graphs
        plt.figure(run)
        plt.plot(c_p[run], y_values_mm)
        plt.suptitle("Coefficient of Pressure vs y Distance from Chord Line")
        plt.title(f"AOA: {-4 + run*4} degrees")
        plt.ylabel("y distance (mm)")
        plt.xlabel("c_p")
        plt.grid()
        plt.show()

    plt.figure(5)
    plt.plot(np.rad2deg(angles_of_attack), c_d, label="Lab 6: Pressure Wake Calculation")
    plt.plot(np.rad2deg(angles_of_attack), lab_5_c_d_data, label="Lab 5: Surface Pressure Calculation")
    plt.suptitle("Coefficient of Drag on Airfoil")
    plt.title("Comparison of Lab 5 and Lab 6 Methods")
    # plt.ylim([-.15, 1.8])
    plt.xlabel("AOA (degrees)")
    plt.ylabel("c_d")
    plt.legend()
    plt.grid()
    plt.show()


################################################################################################################
# Hotwire Anemometer Calibration

PITOT_TUBE_CALIB = 1/0.0016
HOTWIRE_DATA_PATH = os.path.join(os.getcwd(), 'AER E 344/Lab 6/HotWire Data')

# Stores data from hotwire calibration process
hotwire_data = pd.read_csv(os.path.join(HOTWIRE_DATA_PATH, "HotWire Student Data.csv"), encoding="utf-8")

# Stores voltages measured from pitot tube
pitot_tube_voltages = hotwire_data.values[:,1]
# Stores velocities computed from pitot tube voltages
pitot_tube_velocities = []
# Computes velocities per voltage percentage
for run in range(len(pitot_tube_voltages)):
    temp_pressure = PITOT_TUBE_CALIB * (pitot_tube_voltages[run] - pitot_tube_voltages[0])
    pitot_tube_velocities.append(np.sqrt(temp_pressure * 2 / DENSITY))

# Stores voltages measured from hotwire    
hotwire_voltages = hotwire_data.values[:,2]

# Computes 4th order fit for hotwire voltage vs airspeed
a, b, c, d, e = np.polyfit(hotwire_voltages, pitot_tube_velocities, 4)
# Function for hot
hotwire_fun = lambda x, a, b, c, d, e : a*x**4 + b*x**3 + c*x**2 + d*x + e
# Span for x
x_span = np.linspace(0.958,1.583,1000)
coefficient_of_determination = r2_score(pitot_tube_velocities, hotwire_fun(hotwire_voltages, a, b, c, d, e))

plt.plot(hotwire_voltages, pitot_tube_velocities, 'o', label="x_i : Recorded Data")
plt.plot(hotwire_voltages, hotwire_fun(hotwire_voltages, a, b, c, d, e), label="p(x_i)")
plt.plot(x_span, hotwire_fun(x_span, a, b, c, d, e), '--', label="p(x)")
plt.suptitle("Hotwire Anemometer Calibration Curve for 0-20 m/s")
plt.title(f'p(x) = {a:.1f}x^4 + {b:.1f}x^3 + {c:.1f}x^2 + {d:.1f}x + {e:.1f}')
plt.text(1, 17, f'R^2 = {coefficient_of_determination:.3f}', fontweight='bold', bbox=dict(facecolor='Green', alpha=0.3))
plt.xlabel("Hotwire Recorded Voltage")
plt.ylabel("Airspeed")
plt.legend(loc="lower right")
plt.grid()
plt.show()
