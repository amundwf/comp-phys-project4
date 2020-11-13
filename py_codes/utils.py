# This file contains useful functions to be used in other .py files.

import numpy as np
import matplotlib.pyplot as plt
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
    data = pd.DataFrame(data, columns=["MC_cycles", "E", "M"])
    
    f = plt.figure(figsize=(18, 10), dpi=80, facecolor='w', edgecolor='k')
    
    plt.subplot(211)
    plt.plot(data['MC_cycles'], data['E']/4, '.')
    plt.xlabel("Monte Carlo cycles")
    plt.ylabel("Expectation Energy")
    plt.title("Expectation values for {} Monte Carlo cycles for a 2x2 spin matrix".format(len(data["MC_cycles"])))
    
    plt.subplot(212)
    plt.plot(data['MC_cycles'], data['M']/4, '.')
    plt.xlabel("Monte Carlo cycles")
    plt.ylabel("Expectation Magnetisation")
    
plot_E_and_M()