// The solver class for running the Ising model.
#ifndef ISING_SOLVER_HPP
#define ISING_SOLVER_HPP

#include "utils.hpp"
#include <iostream>
#include <cmath>
#include <random>
#include <armadillo>

class IsingSolver{
private:
    double L; // This is the number of rows and columns of the spin matrix.
    // Assuming quadratic (LxL) spin matrix.

    double T; // The current temperature of the system
    double kB = 1; // The Boltzmann constant
    //double J;
    int J = 1; // The exchange energy between spin neighbours. We will only use J=1,
    // so it is set to the type 'int'.

    int N_MC; // Chosen number of Monte Carlo cycles to run the metropolis algorithm.

    // The solver needs to calculate these quantities:
    double E;       // The current energy of the system.
    int M;          // The current net magnetization M of the system.
    
    // Lists of energy and magnetisation over the MC cycles:
    arma::vec E_list; // Energy values for each MC cycle.
    arma::vec M_list; // Magnetisation values for each MC cycle.
    arma::vec E2_list; // Values of E^2 for each MC cycle.
    arma::vec M2_list; // Values of M^2 for each MC cycle.
    arma::vec M_abs_list; // Values of |M| for each MC cycle.

    // Quantities based on expected values:
    double E_mean;  // The mean energy of the system over the Monte Carlo cycles.
    double M_mean;  // The mean net magnetization of the system.
    double M_abs_mean; // // The mean net absolute value of the magnetization of the system.
    double C_V;     // The specific heat of the system.
    double chi;     // The susceptibility of the system.
    
    arma::imat spinMatrix;  // This is the current spin lattice.
    arma::imat PBC_spinMatrix; // N+2 x N+2 matrix (containing periodic boundary
    // conditions along the edges of the matrix).

    void make_PBC_spinMatrix(); // Makes the PBC matrix from spinMatrix.

    arma::Col<int> dE_values; // Contains the 5 possible values of dE.
    arma::vec weights; // Contains the corresponding 5 possible values of e^{-beta*dE}.
    void set_dE_values_and_weights();
public:
    // Constructor with a spin state matrix as input:
    IsingSolver(arma::imat spinMatrix, double T, int N_MC);
    // Constructor that creates a random spin array, with chosen matrix dimensionality L
    // as input:
    IsingSolver(int L, double T, int N_MC);

    double calculate_E(); // Calculates the current energy of the lattice (PBC).
    int calculate_M(); // Calculates the current net magnetization |M| of the lattice.

    void calculate_mean_results(); // Calculates E_mean, M_mean, C_V from E_list and M_list.

    arma::imat get_spinMatrix(); // Output the spin matrix
    void print_spinMatrix();
    void print_PBCSpinMatrix();

    void print_E_list_and_M_list();
    arma::mat get_E_list_M_list();

    arma::Row<double> get_mean_results();

    // Running the metropolis algorithm:
    void metropolis_one_time(); // Run the Metropolis algorithm one time (one MC cycle?)
    void run_metropolis_full(); // Run the Metropolis algorithm for N_MC MC cycles
};
#endif