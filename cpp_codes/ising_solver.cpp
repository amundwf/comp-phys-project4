// The solver class for running the Ising model.

#include "ising_solver.hpp"
#include <random>

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

// Constructor with an initial spin matrix as input:
IsingSolver::IsingSolver(imat spinMatrix, double J){
    this->spinMatrix = spinMatrix;
	this->J = J;
    L = spinMatrix.n_rows;

    // Set the PBC matrix (by using the newly set spinMatrix):
    make_PBC_spinMatrix();

    // Calculate the energy for this spin matrix:
    //E = calculate_E();
    // Calculate the magnetization for this spin matrix:
	//M = calculate_M();
}

// Constructor which initiates a random spin matrix:
IsingSolver::IsingSolver(int n, double J){
    this->J = J;
	L = n;
	mat randMatrix = randu(n,n); // Uniform random elements between 0 and 1.
    
    double elem;
    // Make the matrix binary (-1 and 1):
    for (int i=0; i<=n-1; i++){
        for (int j=0; j<=n-1; j++){
            elem = randMatrix(i,j);
            // Elements < 0.5 are set to -1 and elements > 0.5 are set to +1.
            if (elem < 0.5) randMatrix(i,j) = -1;
            else randMatrix(i,j) = 1;
        }
    }
    // Convert to integer matrix:
    spinMatrix = conv_to<imat>::from(randMatrix);

    // Set the PBC matrix (by using the newly set spinMatrix):
    make_PBC_spinMatrix();
}


//double IsingSolver::calculate_E(){
    // This function calculates the current energy of the spin system.


/*
  double energy = 0.0;
  for (int i = 0; i < L-1; i++){
    for (int j = 0; j < L; j++){
      //Sums through the vertical contributions of the spin interaction
      energy += lat[i][j] * lat[i+1][j];
    }
  }

  for (int i = 0; i < L; i++){
    for (int j = 0; j < L-1; j++){
      //Sums through the horizontal contributions of the spin interactions
      energy += lat[i][j] * lat[i][j+1];
    }
  }

  for (int i = 0; i < L; i++){
    //PBC vertical edges
    energy += lat[i][0]*lat[i][L-1];
  }

  for (int i = 0; i < L; i++){
    //PBC horizontal edges
    energy += lat[0][i]*lat[L-1][i];
  }

  return -energy;
}
*/

void IsingSolver::print_spinMatrix(){
    // This function prints the current spin system (the spinSystem
    // matrix).
    spinMatrix.print("spinMatrix:");
}

void IsingSolver::print_PBCSpinMatrix(){
    // This function prints the current PBC_spinMatrix.
    PBC_spinMatrix.print("PBC_spinMatrix:");
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