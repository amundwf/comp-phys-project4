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

    L = spinMatrix.n_rows;

    // Set the PBC matrix (by using the newly set spinMatrix):
    make_PBC_spinMatrix();
    // Set the 5 dE values and corresponding weights (Boltzmann PDE):
    set_dE_values_and_weights();

    // Calculate the energy for this spin matrix:
    E = calculate_E();
    // Calculate the magnetization for this spin matrix:
	M = calculate_M();

    // Instantiate E_list and M_list:
    E_list = vec(N_MC+1, fill::zeros); M_list = Col<int>(N_MC+1, fill::zeros);
    // Add E and M to the first elements in their lists (as function of MC cycles):
    E_list(0) = E; M_list(0) = M;
}

// Constructor which initiates a random spin matrix:
IsingSolver::IsingSolver(int L, double T, int N_MC){
    // Assign the chosen values to the IsingSolver object:
    this->L = L;
    this->T = T;
    this->N_MC = N_MC;

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
    E_list = vec(N_MC+1, fill::zeros); M_list = Col<int>(N_MC+1, fill::zeros);
    // Add E and M to the first elements in their lists (as function of MC cycles):
    E_list(0) = E; M_list(0) = M;
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

void IsingSolver::run_metropolis_full(){
    // This function runs metropolis_one_time() many times (N_MC number of times).

    for (int i=1; i<=N_MC; i++){ // N_MC times:
        // Run one MC cycle:
        metropolis_one_time();
        
        // Add the new energy and magnetisation to the lists of E and M (as
        // functions of # MC cycles):
        E_list(i) = E; M_list(i) = M;
    }
}