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

runCppCode = False
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
data = pd.DataFrame(data, columns=["MC_cycle", "E", "M", "E2", "M2", "M_abs"])
MC_cycle_list = data["MC_cycle"]
E_list = data["E"]
M_list = data["M"]
E2_list = data["E2"]
M2_list = data["M2"]
M_abs_list = data["M_abs"]

# Mean value quantities:
filename = "4c_MVs.csv"
filePath = os.path.join(directory, filename)
data = np.loadtxt(filePath, skiprows=1, delimiter=",")
data = pd.DataFrame([data], columns=["E_mean", "M_mean", "M_abs_mean", "C_V", "chi"])
E_mean = (data["E_mean"])[0]
M_mean = (data["M_mean"])[0]
M_abs_mean = (data["M_abs_mean"])[0]
C_V = (data["C_V"])[0]
chi = (data["chi"])[0]

print(); print("Mean value quantities:")
print("E_mean: ", E_mean)
print("M_mean: ", M_mean)
print("M_abs_mean: ", M_abs_mean)
print("C_V: ", C_V)
print("chi: ", chi)



# if plotCodeWord == "E_list":
# (plot M_list)
plt.plot(MC_cycle_list, E_list)


# else if plotCodeWord == "M_list":
# (plot E_list)


plt.show()