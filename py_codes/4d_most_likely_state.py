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

runCppCode = 0
# Set to false if you have the results files and just want to plot the results.
if runCppCode == True: 
    # Compile and run the C++ files (this is exactly what is in the makefile):
    os.system("echo compiling C++ codes...")
    #os.system("g++ -o main.out ../cpp_codes/main.cpp ../cpp_codes/utils.cpp ../cpp_codes/planet.cpp ../cpp_codes/solver.cpp -larmadillo")
    os.system("g++ -O3 -o main.out ../cpp_codes/main.cpp ../cpp_codes/utils.cpp ../cpp_codes/ising_solver.cpp -larmadillo")
    os.system("echo executing...")
    os.system("./main.out")

# The results directory:
directory = "../results/4d/"

# Mean values as function of # MC cycles:
filename = "4d_mean_values_vs_MC_cycles.csv"
filePath = os.path.join(directory, filename) # The full file path.
data = np.loadtxt(filePath, skiprows=1, delimiter=",")
# Get the columns of data as lists:
data = pd.DataFrame(data, columns=["MC_cycle", "E_mean", "E2_mean", "M_mean", \
    "M_abs_mean", "M2_mean", "C_V", "chi"])
MC_cycle_list = data["MC_cycle"].tolist() # tolist() converts dataframe to a list.
E_mean_list = data["E_mean"].tolist()
E2_mean_list = data["E2_mean"].tolist()
M_mean_list = data["M_mean"].tolist()
M_abs_mean_list = data["M_abs_mean"].tolist()
M2_mean_list = data["M2_mean"].tolist()
C_V_list = data["C_V"].tolist()
chi_list = data["chi"].tolist()

# E_list and M_list:
filename = "4d_E_list_M_list.csv"
filePath = os.path.join(directory, filename) # The full file path.
data = np.loadtxt(filePath, skiprows=1, delimiter=",")
# Get the columns of data as lists:
data = pd.DataFrame(data, columns=["MC_cycle", "E", "E2", "M", "M_abs", "M2"])
MC_cycle_list2 = data["MC_cycle"].tolist() # tolist() converts dataframe to a list.
E_list = data["E"].tolist()
E2_list = data["E2"].tolist()
M_list = data["M"].tolist()
M_abs_list = data["M_abs"].tolist()
M2_list = data["M2"].tolist()


# Plot the results:
# Mean values:
plotCodeWord = "E_mean"
#plotCodeWord = "E2_mean"
#plotCodeWord = "M_mean"
#plotCodeWord = "M_abs_mean"
#plotCodeWord = "M2_mean"
#plotCodeWord = "C_V"
#plotCodeWord = "chi"
# Or plot the lists:
#plotCodeWord = "E_list"
#plotCodeWord = "M_list"

labelSize = 13
titleSize = 12
markerSize = 5

if plotCodeWord == "E_mean":
    plt.plot(MC_cycle_list, E_mean_list, 'r.', label='E_mean', markersize=markerSize)
    #plt.ylabel(r'$\langleE\rangle$', fontsize=labelSize)
    plt.ylabel('E_mean', fontsize=labelSize)
    plt.suptitle('Mean energy of the 20x20 spin system vs. Monte Carlo cycles', fontsize=titleSize)

elif plotCodeWord == "E2_mean":
    plt.plot(MC_cycle_list, E2_mean_list, 'r.', label='E2_mean', markersize=markerSize)
    plt.ylabel('E2_mean', fontsize=labelSize)
    plt.suptitle('Mean of ' + r'$E^2$' + ' for the 20x20 spin system vs. Monte Carlo cycles', fontsize=titleSize)

elif plotCodeWord == "M_mean":
    plt.plot(MC_cycle_list, M_mean_list, 'r.', label='M_mean', markersize=markerSize)
    plt.ylabel('M_mean', fontsize=labelSize)
    plt.suptitle('Mean of ' + r'$M$' + ' for the 20x20 spin system vs. Monte Carlo cycles', fontsize=titleSize)

elif plotCodeWord == "M_abs_mean":
    plt.plot(MC_cycle_list, M_mean_list, 'r.', label='|M|_mean', markersize=markerSize)
    plt.ylabel('|M|_mean', fontsize=labelSize)
    plt.suptitle('Mean of ' + r'$|M|$' + ' for the 20x20 spin system vs. Monte Carlo cycles', fontsize=titleSize)

elif plotCodeWord == "M2_mean":
    plt.plot(MC_cycle_list, M2_mean_list, 'r.', label='M2_mean', markersize=markerSize)
    plt.ylabel('M2_mean', fontsize=labelSize)
    plt.suptitle('Mean of ' + r'$M^2$' + ' for the 20x20 spin system vs. Monte Carlo cycles', fontsize=titleSize)

elif plotCodeWord == "C_V":
    plt.plot(MC_cycle_list, C_V_list, 'r.', label=r'$C_V$', markersize=markerSize)
    plt.ylabel(r'$C_V$', fontsize=labelSize)
    plt.suptitle('Numerical heat capacity for the 20x20 spin system vs. Monte Carlo cycles', fontsize=titleSize)

elif plotCodeWord == "chi":
    plt.plot(MC_cycle_list, chi_list, 'r.', label=r'$\chi$', markersize=markerSize)
    plt.ylabel(r'$\chi$', fontsize=labelSize)
    plt.suptitle('Numerical magnetic susceptibility for the 20x20 spin system vs. Monte Carlo cycles', fontsize=titleSize)


elif plotCodeWord == "E_list":
    # Plot E_list
    plt.plot(MC_cycle_list2, E_list, '.', label=r'$E$', markersize=markerSize)
    # Plot the final mean value as well (as a straight horizontal line):
    plt.plot([MC_cycle_list2[0], MC_cycle_list2[-1]], [E_mean_list[-1], E_mean_list[-1]], 'r-', label='E_mean')
    plt.ylabel(r'$E$', fontsize=labelSize)
    plt.suptitle('The energy of the 20x20 spin system vs. Monte Carlo cycles', fontsize=titleSize)

elif plotCodeWord == "M_list":
    # Plot M_list
    plt.plot(MC_cycle_list2, M_list, '.', label=r'$M$', markersize=markerSize)
    # Plot the final mean value as well:
    plt.plot([MC_cycle_list2[0], MC_cycle_list2[-1]], [M_mean_list[-1], M_mean_list[-1]], 'r-', label='M_mean')
    plt.ylabel(r'$M$', fontsize=labelSize)
    plt.suptitle('The magnetisation of the 20x20 spin system vs. Monte Carlo cycles', fontsize=titleSize)


plt.xlabel(r'Number of Monte Carlo cycles', fontsize=labelSize)
plt.legend()
plt.grid()
plt.show()