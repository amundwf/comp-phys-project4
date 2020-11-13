# This file contains useful functions to be used in other .py files.

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})
import pandas as pd

def energy_EV_analytical_2x2(T, J, kB):
    # Returns the analytical expectation value for the 2x2 case,
    # as a function of temperature and the interaction term J (and
    # the Boltzmann constant kB).
    beta = 1/(kB*T)
    x = 8*J*beta
    print("x: "); print(x)
    print(np.sinh(x)); print(np.cosh(x))
    # The analytical expression:
    return -8*J*(np.sinh(x) / (np.cosh(x) + 3))

'''
def read_matrix_from_CSV(filename, directory):
    directory = "../results/3g_three_body/"

    # Read the planet names from the planet_names.csv file:
    
    #planetNamesData = np.loadtxt(directory + "planet_names.csv", skiprows=0, delimiter=",")
    #planetNamesData = pd.DataFrame(planetNamesData, columns=["planet_names"])
    #planetNames = planetNamesData["planet_names"]
    
    # Read in the file as a dataframe:
    dataFrame = pd.read_csv(directory + filename, header=None)
'''

def plot_spinMatrix(spinMatrix):
    # input: spinMatrix: A numpy array containing values -1 and 1.
    # This function plots a colormap of black (-1) and white (+1) squares.

    tempColorMatrix = spinMatrix
    tempColorMatrix[spinMatrix==-1] = 0
    plt.imshow(tempColorMatrix, cmap = 'gray') # The 'gray' colormap is binary:
    # black for values below 1 and white for values >= 1.
    plt.show()

def plot_E_and_M():
    
    directory = "../results/4c_ising/"
    filename = "blabla_4c.csv"
    
    data = np.loadtxt(directory +filename, skiprows=1, delimiter=",")
    data = pd.DataFrame(data, columns=["MC_cycles", "E", "E2", "M", "M2", "Mabs"])
    
    T = 1.0
    L2 = 2*2
    MC = data["MC_cycles"]
    N_MC = len(data['MC_cycles']) -1
    E = data["E"]
    E2 = data["E2"]
    M = data["M"]
    M2 = data["M2"]
    Mabs = data["Mabs"]
    
    CV = (E2 - E*E)/T*T
    chi = (M2 - M*M)/T
    
    f = plt.figure(figsize=(18, 10), dpi=80, facecolor='w', edgecolor='k')
    
    plt.subplot(221)
    plt.plot(MC, E/L2, '.')
    plt.xscale("log")
    plt.xlabel("Monte Carlo cycles")
    plt.ylabel("Expectation Energy")
    
    plt.subplot(222)
    plt.plot(MC, Mabs/L2, '.')
    plt.xscale("log")
    plt.xlabel("Monte Carlo cycles")
    plt.ylabel("Expectation Magnetisation")
    
    plt.subplot(223)
    plt.plot(MC, CV , '.')
    plt.xscale("log")
    plt.xlabel("Monte Carlo cycles")
    plt.ylabel("Specific Heat")
    
    plt.subplot(224)
    plt.plot(MC, chi , '.')
    plt.xscale("log")
    plt.xlabel("Monte Carlo cycles")
    plt.ylabel("Magnetic Susceptibility")
    
    f.suptitle("Expectation values for {} Monte Carlo cycles for a 2x2 spin matrix".format(N_MC))
    
    
    
plot_E_and_M()