# This script runs the Metropolis algorithm for the 2x2 spin system.
# It calculates and plots the energy and magnetisation as a function of
# Monte Carlo cycles.
# It also calculates the mean energy, mean magnetisation, heat capacity
# and susceptibility.
# Optional: This program can also plot the analytical values for the
# 2x2 case.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import os
import utils as ut

runCppCode = 1
# Set to false if you have the results files and just want to plot the results.
if runCppCode == True: 
    # Compile and run the C++ files (this is exactly what is in the makefile):
    os.system("echo compiling C++ codes...")
    #os.system("g++ -o main.out ../cpp_codes/main.cpp ../cpp_codes/utils.cpp ../cpp_codes/planet.cpp ../cpp_codes/solver.cpp -larmadillo")
    os.system("g++ -O3 -o main.out ../cpp_codes/main.cpp ../cpp_codes/utils.cpp ../cpp_codes/ising_solver.cpp -larmadillo")
    os.system("echo executing...")
    os.system("./main.out")


# The results directory:
directory = "../results/4c_ising/"

# Read the planet names from the planet_names.csv file:
'''
planetNamesData = np.loadtxt(directory + "planet_names.csv", skiprows=0, delimiter=",")
planetNamesData = pd.DataFrame(planetNamesData, columns=["planet_names"])
planetNames = planetNamesData["planet_names"]
'''

# E_list and M_list:
filename = "4c_E_list_M_list.csv"
filePath = os.path.join(directory, filename) # The full file path.
data = np.loadtxt(filePath, skiprows=1, delimiter=",")
# Get the columns of data as lists:
data = pd.DataFrame(data, columns=["MC_cycle", "E", "E2", "M", "M_abs", "M2"])
MC_cycle_list = data["MC_cycle"].tolist() # tolist() converts dataframe to a list.
E_list = data["E"].tolist()
E2_list = data["E2"].tolist()
M_list = data["M"].tolist()
M_abs_list = data["M_abs"].tolist()
M2_list = data["M2"].tolist()


# Mean value quantities:
filename = "4c_MVs.csv"
filePath = os.path.join(directory, filename)
data = np.loadtxt(filePath, skiprows=1, delimiter=",")
data = pd.DataFrame([data], columns=["E_mean", "E2_mean", "M_mean", "M_abs_mean", "M2_mean", "C_V", "chi"])
E_mean = (data["E_mean"].tolist())[0]
E2_mean = (data["E2_mean"].tolist())[0]
M_mean = (data["M_mean"].tolist())[0]
M_abs_mean = (data["M_abs_mean"].tolist())[0]
M2_mean = (data["M2_mean"].tolist())[0]
C_V = (data["C_V"].tolist())[0]
chi = (data["chi"].tolist())[0]

print(); print("Numerical mean value quantities:")
print("E_mean: ", E_mean)
print("E2_mean: ", E2_mean)
print("M_mean: ", M_mean)
print("M_abs_mean: ", M_abs_mean)
print("M2_mean: ", M2_mean)
print("C_V: ", C_V)
print("chi: ", chi)

print(); print("Analytical mean value quantities:")
J = 1.0; kB = 1.0; T = 1.0
E_mean_analytical = ut.E_mean_analytical_2x2(J, kB, T)
E2_mean_analytical = ut.E2_mean_analytical_2x2(J, kB, T)
M_abs_mean_analytical = ut.M_abs_mean_analytical_2x2(J, kB, T)
M2_mean_analytical = ut.M2_mean_analytical_2x2(J, kB, T)
C_V_analytical = ut.C_V_analytical_2x2(J, kB, T)
chi_analytical = ut.chi_analytical_2x2(J, kB, T)
print("E_mean: ", E_mean_analytical)
print("E2_mean: ", E2_mean_analytical)
print("M_mean: 0?")
print("M_abs_mean: ", M_abs_mean_analytical)
print("M2_mean: ", M2_mean_analytical)
print("C_V: ", C_V_analytical)
print("chi: ", chi_analytical)


#plotCodeWord = "E_list"
plotCodeWord = "M_list"

labelSize = 13
titleSize = 12

if plotCodeWord == "E_list":
    # Plot E_list
    plt.plot(MC_cycle_list, E_list, '.', label='E_list', markersize=3)
    # Plot the mean as well (as a straight horizontal line)
    plt.plot([MC_cycle_list[0], MC_cycle_list[-1]], [E_mean, E_mean], 'r-', label='E_mean')
    #plt.plot([MC_cycle_list[0], MC_cycle_list[-1]], [E_mean, E_mean])
    plt.ylabel(r'$E$', fontsize=labelSize)
    plt.suptitle('Energy of 2x2 spin system vs. Monte Carlo cycles', fontsize=titleSize)

elif plotCodeWord == "M_list":
    # Plot M_list
    plt.plot(MC_cycle_list, M_list, '.')
    plt.ylabel(r'$M$', fontsize=labelSize)
    plt.suptitle('Magnetisation of 2x2 spin system vs. Monte Carlo cycles', fontsize=titleSize)

plt.xlabel(r'Number of Monte Carlo cycles', fontsize=labelSize)

#plt.xlim(-1.5, 1.5)
#plt.ylim(-1.5, 1.5)
#plt.axis('equal')
plt.legend()
plt.grid()
plt.show()