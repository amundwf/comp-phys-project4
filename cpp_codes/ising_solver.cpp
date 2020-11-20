// The solver class for running the Ising model.

#include "ising_solver.hpp"
#include <random>
#include <time.h>


using namespace std;
using namespace arma;

void IsingSolver::make_PBC_spinMatrix(){
    // This function is used inside the constructor IsingSolver(imat spinMatrix).
    // In that constructor, squareMatrix and L are set before this function is called,
    // so they have correct values.

    // Construct the PBC spin matrix:
    PBC_spinMatrix = imat(L+2, L+2, fill::zeros);

    // Copy the square matrix into the center:
    PBC_spinMatrix(span(1,L), span(1,L)) = spinMatrix;
    // Fill each of the 4 edge sides of the PBC matrix with the opposite end
    // side of the square matrix:
    // Left vertical edge side:
    PBC_spinMatrix(span(1,L), 0) = spinMatrix(span(0,L-1), L-1);
    // Right vertical edge side:
    PBC_spinMatrix(span(1,L), L+1) = spinMatrix(span(0,L-1), 0);
    // Upper horizontal edge side:
    PBC_spinMatrix(0, span(1,L)) = spinMatrix(L-1, span(0,L-1));
    // Lower horizontal edge side:
    PBC_spinMatrix(L+1, span(1,L)) = spinMatrix(0, span(0,L-1));
}

void IsingSolver::set_dE_values_and_weights(){
    // This function sets dE_values and pre-calculates weights.
    // It is called in IsingSolver constructors.
    dE_values = Col<int>("-8, -4, 0, 4, 8");

    weights = vec(17); //8-(-8)+1 = 17 elements.
    // Pre-calculate w:
    double beta = 1/(kB*T);
    for (int idx=0; idx<=16; idx++){
        int dE = idx-8; // Transforms index to corresponding dE
        weights(idx) = exp(-beta*dE);
    }
    // Usage: When performing a spin flip and calculating dE caused by that
    // flip, the corresponding weight can be accessed by transforming dE to
    // an index of weights by adding 8 to dE:
    // int dE = ...; w = weights(dE+8); 
}

// Constructor with an initial spin matrix as input:
IsingSolver::IsingSolver(imat spinMatrix, double T, int N_MC){
    // Assign the chosen values to the IsingSolver object:
    this->spinMatrix = spinMatrix;
    this->T = T;
    this->N_MC = N_MC;
    // Equilibrium cycles subtracted from full MC cycles.
    this->N1_MC = 0.9*N_MC;
    
    L = spinMatrix.n_rows;
    this->L2 = L*L;
    // Set the PBC matrix (by using the newly set spinMatrix):
    make_PBC_spinMatrix();
    // Set the 5 dE values and corresponding weights (Boltzmann PDE):
    set_dE_values_and_weights();

    // Calculate the energy for this spin matrix:
    E = calculate_E();
    // Calculate the magnetization for this spin matrix:
	M = calculate_M();

    // Instantiate E_list and M_list:
    E_list = vec(N_MC+1, fill::zeros);
    M_list = vec(N_MC+1, fill::zeros);
    M_abs_list = vec(N_MC+1, fill::zeros);

    // Initialise E2_list and M2_list:
    E2_list = vec(N_MC+1, fill::zeros);
    M2_list = vec(N_MC+1, fill::zeros);
    
    // Add E and M to the first elements in their lists (as function of MC cycles):
    E_list(0) = E; M_list(0) = M;
    E2_list(0) = E*E; M2_list(0) = M*M;
    M_abs_list(0) = fabs(M);

    // Number of accepted flips per MC cycle.
    flipsAccepted_list = vec(N_MC+1, fill::zeros);
}

// Constructor which initiates a random spin matrix:
IsingSolver::IsingSolver(int L, double T, int N_MC){
    // Assign the chosen values to the IsingSolver object:
    this->L = L;
    this->T = T;
    this->N_MC = N_MC;
    // Equilibrium cycles subtracted from full MC cycles.
    this->N1_MC = 0.9*N_MC;
    this->L2 = L*L;

    arma_rng::set_seed_random(); // Set seed to a random value.
	mat randMatrix = randu(L,L); // Uniform random elements between 0 and 1.

    // Make the matrix binary (-1 and 1):
    for (int i=0; i<=L-1; i++){
        for (int j=0; j<=L-1; j++){
            double elem = randMatrix(i,j);
            // Elements < 0.5 are set to -1 and elements > 0.5 are set to +1.
            if (elem < 0.5) randMatrix(i,j) = -1;
            else randMatrix(i,j) = 1;
        }
    }
    // Convert to integer matrix:
    spinMatrix = conv_to<imat>::from(randMatrix);

    // Set the PBC matrix (by using the newly set spinMatrix):
    make_PBC_spinMatrix();
    // Set the 5 dE values and corresponding weights (Boltzmann PDE):
    set_dE_values_and_weights();

    // Calculate the energy for this spin matrix:
    E = calculate_E();
    // Calculate the magnetization for this spin matrix:
	M = calculate_M();

    // Instantiate E_list and M_list:
    E_list = vec(N_MC+1, fill::zeros);
    M_list = vec(N_MC+1, fill::zeros);
    M_abs_list = vec(N_MC+1, fill::zeros);

    // Initialise E2_list and M2_list:
    E2_list = vec(N_MC+1, fill::zeros);
    M2_list = vec(N_MC+1, fill::zeros);
    
    // Add E and M to the first elements in their lists (as function of MC cycles):
    E_list(0) = E; M_list(0) = M;
    E2_list(0) = E*E; M2_list(0) = M*M;
    M_abs_list(0) = fabs(M);

    flipsAccepted_list = vec(N_MC+1, fill::zeros);

    E_total = 0.0;
    E2_total = 0.0;
    M_total = 0.0;
    M_abs_total = 0.0;
    M2_total = 0.0;
}

imat IsingSolver::get_spinMatrix(){
    return spinMatrix;
}

void IsingSolver::print_spinMatrix(){
    // This function prints the current spin system (the spinSystem
    // matrix).
    spinMatrix.print("spinMatrix:");
}

void IsingSolver::print_PBCSpinMatrix(){
    // This function prints the current PBC_spinMatrix.
    PBC_spinMatrix.print("PBC_spinMatrix:");
}

void IsingSolver::print_E_list_and_M_list(){
    // This function prints E_list and M_list in an array with two columns.

    vec M_list_double = conv_to<vec>::from(M_list); // Convert M_list to double vector
    // Paste E_list and M_list together in an array with two columns:
    mat E_M_array = mat(N_MC+1, 2, fill::zeros);
    //E_list.print();
    //M_list_double.print();
    
    E_M_array(span(0,N_MC), 0) = E_list; // Energy in first column
    E_M_array(span(0,N_MC), 1) = M_list_double; // Magnetisation in second column
    // Print the array:
    E_M_array.print("E_list and M_list:");
}

mat IsingSolver::get_E_list_M_list(){
    // This function returns the results ,should rename it. Also 
    // it divides by the number of spins. 
    vec MC_list = vec(N_MC+1);
    for (int i=0; i<=N_MC; i++){
        MC_list(i) = i;
    }

    vec M_list_double = conv_to<vec>::from(M_list); // Convert M_list to double vector
    // Paste E_list and M_list together in an array with two columns:
    mat E_M_array = mat(N_MC+1, 7, fill::zeros);

    E_M_array(span(0,N_MC), 0) = MC_list;
    E_M_array(span(0,N_MC), 1) = E_list;
    E_M_array(span(0,N_MC), 2) = E2_list;
    E_M_array(span(0,N_MC), 3) = M_list_double;
    E_M_array(span(0,N_MC), 4) = M_abs_list;
    E_M_array(span(0,N_MC), 5) = M2_list;
    E_M_array(span(0,N_MC), 6) = flipsAccepted_list;

    // Return the array:
    return E_M_array;
}

int IsingSolver::calculate_M(){
	int M = 0; // The magnetization.
	for (int i=0; i<=L-1; i++){
		for (int j=0; j<=L-1; j++){
			M += spinMatrix(i,j);
		}
	}
	//M = abs(M);
	this->M = M;
	return M;
}

double IsingSolver::calculate_E(){
	// This function calculates the energy of the current spin lattice, replaces
	// the energy member variable E, and also returns this energy.
	// Go through each element in the PBC spin matrix. 
	double E = 0;
	double E_currentSpin = 0; // Temp variable for the energy of each spin.

	// Remember that the PBC spin matrix has slightly different indices compared
	// to spinMatrix due to its increased size ((L+2)x(L+2)).
	int s1; // The spin value of the current spin in a row or column.
	int s0; // The spin value of the previous spin (in a row or a column).

	// The total energy will be calculated by first getting the energy from the row
	// strings of spins and then getting the energy from the column strings of spins:

	// Get the energy sum of all of the row strings:
	for (int r=1; r<=L; r++){ // For each row, get the binding energies between the
	// spins in that row.
		for (int c=1; c<=L; c++){ // Column number. Include the last element (index
		// L+1) for the periodic boundary conditions.
			s1 = PBC_spinMatrix(r,c);
			s0 = PBC_spinMatrix(r,c-1); // Previous spin in row
			E += s1*s0;
		}
	}
	// Get the energy sum of all of the column strings:
	for (int c=1; c<=L; c++){ // For each column, get the binding energies between the
	// spins in that column.
		for (int r=1; r<=L; r++){ // Row number
			s1 = PBC_spinMatrix(r,c);
			s0 = PBC_spinMatrix(r-1,c); // Previous spin in column
			E += s1*s0;
		}
	}

	E *= -J; // Multiply with the factor -J to get the correct energy. 
	this->E = E;
	return E;
}

void IsingSolver::calculate_mean_results(){
    // Calculates E_mean, M_mean, C_V from E_list and M_list and
    // updates these member variables.
    
    // Discard the first 10% of the MC cycles (to reach the most likely state
    // before including the values of E and M in their mean values):
    int N_start = round(double(N_MC)/10);
    vec E_list_shortened = E_list(span(N_start, N_MC));
    vec M_list_shortened = M_list(span(N_start, N_MC));
    vec E2_list_shortened = E2_list(span(N_start, N_MC));
    vec M2_list_shortened = M2_list(span(N_start, N_MC));
    vec M_abs_list_shortened = M_abs_list(span(N_start, N_MC));

    // Now calculate the mean values of these lists:
    E_mean = mean(E_list_shortened);
    M_mean = mean(M_list_shortened);
    E2_mean = mean(E2_list_shortened);
    M2_mean = mean(M2_list_shortened);
    M_abs_mean = mean(M_abs_list_shortened);
    // Heat capacity:
    C_V = (1/(kB*T*T))*(E2_mean - E_mean*E_mean);
    // Susceptibility:
    chi = (1/(kB*T))*(M2_mean - M_abs_mean*M_abs_mean);
}

void IsingSolver::calculate_mean_results_current_MC_cycle(int current_MC_cycle){
    // Same as calculate_mean_results(), but now it includes all MC cycles,
    // that is, it starts at N_start = 1, and only goes up to the current
    // MC cycle. It does not discard the first 10% or so of the cycles.
    // This is for testing purposes, to see when the spin system reaches
    // an equilibrium.
    // int current_MC_cycle: The number of the current MC cycle (1,2,...,N_MC).
    int N_start = 0;
    vec E_list_shortened = E_list(span(N_start, current_MC_cycle));
    vec M_list_shortened = M_list(span(N_start, current_MC_cycle));
    vec E2_list_shortened = E2_list(span(N_start, current_MC_cycle));
    vec M2_list_shortened = M2_list(span(N_start, current_MC_cycle));
    vec M_abs_list_shortened = M_abs_list(span(N_start, current_MC_cycle));

    // Now calculate the mean values of these lists and update the member
    // variables:
    E_mean = mean(E_list_shortened);
    M_mean = mean(M_list_shortened);
    E2_mean = mean(E2_list_shortened);
    M2_mean = mean(M2_list_shortened);
    M_abs_mean = mean(M_abs_list_shortened);
    // Heat capacity:
    C_V = (1/(kB*T*T))*(E2_mean - E_mean*E_mean);
    // Susceptibility:
    chi = (1/(kB*T))*(M2_mean - M_abs_mean*M_abs_mean);
}

Row<double> IsingSolver::get_mean_results(){
    // Returns all four quantities that are based on mean values: <E>, <M>,
    // C_V and chi.
    Row<double> MVs_array(7); // MV = mean value
    MVs_array(0) = E_mean; MVs_array(1) = E2_mean;
    MVs_array(2) = M_mean; MVs_array(3) = M_abs_mean; MVs_array(4) = M2_mean;
    MVs_array(5) = C_V; MVs_array(6) = chi;

    return MVs_array;
}

void IsingSolver::metropolis_one_time(){
    // This function runs the Metropolis algorithm one time (one MC cycle?)
    // (See the metropolis algorithm implementation in the lecture notes p. 438)

    for (int y=0; y<=L-1; y++){ // Go through all spins (but pick a random position
    // each time).
        for (int x=0; x<=L-1; x++){
            // Get a random position in the (PBC) lattice:
            int iy = randrange_int(1, L);
            int ix = randrange_int(1, L);

            // Calculate the energy difference that would be caused by flipping
            // the spin at matrix indices (iy,ix):
            int dE = 2*PBC_spinMatrix(iy,ix)*
                (PBC_spinMatrix(iy,ix-1) + PBC_spinMatrix(iy,ix+1)
                + PBC_spinMatrix(iy-1,ix) + PBC_spinMatrix(iy+1,ix));
                // ^ A sum over the four nearest neighbours.

            // Now perform the Metropolis test using dE:
            double w = weights(dE+8); // Get the weight corresponding to dE (without
            // needing to calculate the exponential).
            double randNum = randrange_float(0,1); // Uniformly random double between 0 and 1.
            if (randNum <= w){ // If the weight is larger than the random number, the spin flip
                // is acceptable. Therefore, perform the spin flip.
                PBC_spinMatrix(iy,ix) *= -1; spinMatrix(iy-1,ix-1) *= -1; 
                // Also update the energy and magnetisation:
                E += dE;
                M += 2*PBC_spinMatrix(iy,ix);
                
            }
        }
    }
}

void IsingSolver::metropolis_one_time_and_update_E_M_lists(int list_idx){
    // Same as metropolis_one_time(), but this function also updates E, E2, M,
    // etc. in E_list, E2_list, M_list, ... respectively, for one index value.
    // int list_idx: Current list index of E_list, E2_list, etc. corresponding to
    // the current MC cycle.
    int flipsAccepted = 0;
    for (int y=0; y<=L-1; y++){ // Go through all spins (but pick a random position
    // each time).
        for (int x=0; x<=L-1; x++){
            // Get a random position in the (PBC) lattice:
            int iy = randrange_int(1, L);
            int ix = randrange_int(1, L);

            // Calculate the energy difference that would be caused by flipping
            // the spin at matrix indices (iy,ix):
            int dE = 2*PBC_spinMatrix(iy,ix)*
                (PBC_spinMatrix(iy,ix-1) + PBC_spinMatrix(iy,ix+1)
                + PBC_spinMatrix(iy-1,ix) + PBC_spinMatrix(iy+1,ix));
                // ^ A sum over the four nearest neighbours.

            // Now perform the Metropolis test using dE:
            double w = weights(dE+8); // Get the weight corresponding to dE (without
            // needing to calculate the exponential).
            double randNum = randrange_float(0,1); // Uniformly random double between 0 and 1.
            if (randNum <= w){ // If the weight is larger than the random number, the spin flip
                // is acceptable. Therefore, perform the spin flip.
                PBC_spinMatrix(iy,ix) *= -1; spinMatrix(iy-1,ix-1) *= -1; 
                // Also update the energy and magnetisation:
                E += dE;
                M += 2*PBC_spinMatrix(iy,ix);
                flipsAccepted += 1;
            }
        }
    }
    // Add the new energy and magnetisation to the lists of E and M (as
    // functions of # MC cycles):
    E_list(list_idx) = E; M_list(list_idx) = M; E2_list(list_idx) = E*E; M2_list(list_idx) = M*M;
    M_abs_list(list_idx)= fabs(M); flipsAccepted_list(list_idx) = flipsAccepted;
}

void IsingSolver::run_metropolis_full(){
    // This function runs metropolis_one_time() many times (N_MC number of times).

    for (int i=1; i<=N_MC; i++){ // N_MC times:
        // Run one MC cycle:
        metropolis_one_time();
        
        // Add the new energy and magnetisation to the lists of E and M (as
        // functions of # MC cycles):
        E_list(i) = E; M_list(i) = M; E2_list(i) = E*E; M2_list(i) = M*M;
        M_abs_list(i)= fabs(M);
    }

    // Now that the metropolis algorithm has been run, calculate
    // and update the expectated value quantities:
    calculate_mean_results();
}

void IsingSolver::init_parallel_variables(){
    // Set all variables to zero. 
    E_total = 0.0;
    E2_total = 0.0;
    M_total = 0.0;
    M_abs_total = 0.0;
    M2_total = 0.0;
}

void IsingSolver::metropolis_one_time_parallel(){
    // Same as metropolis_one_time(), but this function also updates E, E2, M,
    // etc. in E_list, E2_list, M_list, ... respectively, for one index value.
    // int list_idx: Current list index of E_list, E2_list, etc. corresponding to
    // the current MC cycle.
    for (int y=0; y<=L-1; y++){ // Go through all spins (but pick a random position
    // each time).
        for (int x=0; x<=L-1; x++){
            // Get a random position in the (PBC) lattice:
            int iy = randrange_int(1, L);
            int ix = randrange_int(1, L);

            // Calculate the energy difference that would be caused by flipping
            // the spin at matrix indices (iy,ix):
            int dE = 2*PBC_spinMatrix(iy,ix)*
                (PBC_spinMatrix(iy,ix-1) + PBC_spinMatrix(iy,ix+1)
                + PBC_spinMatrix(iy-1,ix) + PBC_spinMatrix(iy+1,ix));
                // ^ A sum over the four nearest neighbours.

            // Now perform the Metropolis test using dE:
            double w = weights(dE+8); // Get the weight corresponding to dE (without
            // needing to calculate the exponential).
            double randNum = randrange_float(0,1); // Uniformly random double between 0 and 1.
            if (randNum <= w){ // If the weight is larger than the random number, the spin flip
                // is acceptable. Therefore, perform the spin flip.
                PBC_spinMatrix(iy,ix) *= -1; spinMatrix(iy-1,ix-1) *= -1; 
                // Also update the energy and magnetisation:
                E += dE;
                M += 2*PBC_spinMatrix(iy,ix);
            }
        }
    }
    // Add the new energy and magnetisation to the lists of E and M (as
    // functions of # MC cycles):
    E_total += E; M_total += M; E2_total += E*E; M2_total += M*M;
    M_abs_total += fabs(M);
}

Row<double> IsingSolver::get_mean_results_parallel(){
    // Returns all four quantities that are based on mean values: <E>, <M>,
    // C_V and chi.

    // The norm will take into account the equil steps. 
    double N_MCdouble = double(N_MC);
    double norm = 1/ (0.9*N_MCdouble);
    double E_mean = E_total*norm;
    double M_mean = M_total*norm;
    double E2_mean = E2_total*norm;
    double M2_mean = M2_total*norm;
    double M_abs_mean = M_abs_total*norm;
    
    double C_V = (1/(kB*T*T))*(E2_mean - E_mean*E_mean);
    // Susceptibility:
    double chi = (1/(kB*T))*(M2_mean - M_abs_mean*M_abs_mean);

    Row<double> MVs_array(7); // MV = mean value
    MVs_array(0) = E_mean/L2; MVs_array(1) = E2_mean/L2;
    MVs_array(2) = M_mean/L2; MVs_array(3) = M_abs_mean/L2; MVs_array(4) = M2_mean/L2;
    MVs_array(5) = C_V/L2; MVs_array(6) = chi/L2;

    return MVs_array;
}