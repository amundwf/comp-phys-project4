# comp-phys-project4
Repository for project 4 in FYS4150 computational physics.

Brief Description: This is a repository for using the 2D Ising model to simulate the ferromagnetic phase transition. 

The script can be run from the make file within cpp_codes; please choose the function you wish to run by uncommenting it in utils.cpp. Or the script can be run from the python files e.g 4c_ising.py within the py_codes folder. The plotting of results can be done directly from the python script e.g 4c_ising.py or from plotting.py.

The results of the simulations are stored in the results folder. 4c_ising is a folder where results are stored for the 2x2 spin matrix, which is compared to the analytic solution in the report. The folder 4d stores results for the expectation values of a 20x20 spin matrix, such as mean energy, magnetisation, specific heat and magnetic susceptibility. Folder 4f stores results for running at several temperatures with many different spin matrices of size 40, 60, 80 and 100. 
