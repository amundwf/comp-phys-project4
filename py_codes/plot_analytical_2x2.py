# Plot the analytical functions for the expectation values of the energy (and more?)

import numpy as np
import matplotlib.pyplot as plt
import utils


#JList = [1e-2, 0.1, 0.5, 1, 2, 4, 5] # Just guessing some random values for J. Plot
# for all of them.
# Use J = 1 and kB = 1 for simplicity. Temperature then has units J. (E = kB*T = 1*T)
JList = [1]

# List of temperature values (in Kelvin) to plot for:
TMin = 0.1; TMax = 2; dT = 0.1; NT = round((TMax-TMin)/dT)
TList = np.linspace(TMin, TMax, NT)

# Boltzmann's constant:
#kB = 1.38064852e-23 # J/K
kB = 1 # Now in units joules.

for i in range(len(JList)):
    J = JList[i]

    #EList = utils.energy_EV_analytical_2x2(TList, J, kB)
    EList = np.zeros(len(TList))
    for j in range(len(TList)):
        T = TList[j]
        E = utils.energy_EV_analytical_2x2(T, J, kB)
        #E_eV = E/1.6e-19
        print(E)
        EList[j] = E

    J_str = format(J, ".3e")
    label = "J = %s" % J_str
    # Plot the energy for this value of J:
    plt.plot(TList, EList, label = label)


labelSize = 13
titleSize = 12
plt.xlabel(r'$T$ (K)', fontsize=labelSize)
plt.ylabel(r'$E$ (eV)', fontsize=labelSize)
plt.legend()
plt.grid()
plt.show()