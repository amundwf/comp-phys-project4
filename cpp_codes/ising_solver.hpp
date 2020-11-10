// The solver class for running the Ising model.
#ifndef ISING_SOLVER_HPP
#define ISING_SOLVER_HPP

#include <iostream>
#include <cmath>
#include <armadillo>

class IsingSolver{
private:
    double T; // The current temperature of the system
    double J; // The exchange energy between spin neighbours

    // The solver needs to calculate these quantities:
    double E;       // The current energy of the system
    double E_mean;  // The mean energy of the system
    int M;          // The current net magnetization M of the system
    double M_mean;  // The mean net magnetization of the system
    double C_v;     // The specific heat of the system
    double Chi;     // The susceptibility of the system

    double L; // This is the number of rows and columns of the spin matrix.
    // Assuming quadratic (LxL) spin matrix.

    arma::imat spinMatrix;  // This is the spins.
    arma::imat PBC_spinMatrix; // N+2 x N+2 matrix (containing periodic boundary
    // conditions along the edges of the matrix).

    int N_MC; // Chosen number of Monte Carlo cycles to run the metropolis algorithm.

    void metropolis_one_time(); // Run the Metropolis algorithm one time (one MC cycle?)
public:
    void make_PBC_spinMatrix();
    // Constructor with a spin state matrix as input:
    IsingSolver(arma::imat spinMatrix, double J);
    // Constructor that creates a random spin array, with matrix dimensionality n (or L)
    // as input:
    IsingSolver(int n, double J);

    double calculate_E(); // Calculates the current energy of the lattice (PBC).
    int calculate_M(); // Calculates the current net magnetization |M| of the lattice.

    void print_spinMatrix();
    void print_PBCSpinMatrix();

};
#endif