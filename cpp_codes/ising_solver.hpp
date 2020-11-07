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

    double L; // This is the number of rows and columns of the spin matrix.
    // Assuming quadratic (LxL) spin matrix.

    arma::imat spinMatrix;  // This is the spins.
    arma::imat PBC_spinMatrix; // N+2 x N+2 matrix (containing periodic boundary
    // conditions along the edges of the matrix).

    int M; // Chosen number of Monte Carlo cycles to run the metropolis algorithm.

    void metropolis_one_time(); // Run the Metropolis algorithm one time (one MC cycle?)
public:
    void make_PBC_spinMatrix();
    // Constructor with a spin state matrix as input:
    IsingSolver(arma::imat spinMatrix);
    // Constructor that creates a random spin array, with matrix dimensionality n
    // as input.
    IsingSolver(int n);

    double calculate_E();
    double calculate_M();

    void printSpinMatrix();
    void printPBCSpinMatrix();

};
#endif