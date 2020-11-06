// The solver class for running the Ising model.
#ifndef ISING_SOLVER_HPP
#define ISING_SOLVER_HPP

#include <iostream>
#include <cmath>
#include <armadillo>

class IsingSolver{
private:
    double T; // Temperature of the system (input variable of the solver?)

    // The solver needs to calculate these:
    //double E_mean; // The mean energy of the system
    //double M_mean; // The mean magnetization |M| of the system
    double E; // The energy of the system
    double M; // The net magnetization |M| of the system
    double C_v; // The specific heat of the system
    double Chi; // The susceptibility of the system

    arma::imat spinMatrix;  // This is the spins.
    arma::imat PBC_spinMatrix; // N+2 x N+2 matrix (containing periodic boundary
    // conditions along the edges of the matrix).
public:
    // Constructor, with a spin state matrix as input:
    IsingSolver(arma::imat spinMatrix);

    void printSpinSystem();

};
#endif