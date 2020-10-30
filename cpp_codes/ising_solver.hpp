// The solver class for running the Ising model.
#ifndef ISING_SOLVER_HPP
#define ISING_SOLVER_HPP

#include <iostream>
#include <cmath>
#include <armadillo>

class IsingSolver{
private:
    double T; // Temperature of the system (input variable)

    // The solver needs to calculate these:
    double E_mean; // The mean energy of the system
    double M_mean; // The mean magnetization |M| of the system
    double C_v; // Specific heat of the system
    double Chi; // Susceptibility of the system
public:

};
#endif