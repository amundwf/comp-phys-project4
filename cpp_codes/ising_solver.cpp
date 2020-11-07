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
IsingSolver::IsingSolver(imat spinMatrix){
    this->spinMatrix = spinMatrix;
    L = spinMatrix.n_rows;

    // Set the PBC matrix (by using the newly set spinMatrix):
    make_PBC_spinMatrix();

    // Calculate the energy for this spin matrix:
    //E = calculate_E();
    // Calculate the magnetization for this spin matrix:
}

// IsingSolver::randomizeSpins(){
    

    // Now update the energy, magnetization, etc.:

//}

// Constructor creating a random spin matrix
IsingSolver::IsingSolver(int n){
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
    L = n;
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

void IsingSolver::printSpinMatrix(){
    // This function prints the current spin system (the spinSystem
    // matrix).
    spinMatrix.print("spinMatrix:");
}

void IsingSolver::printPBCSpinMatrix(){
    // This function prints the current PBC_spinMatrix.
    PBC_spinMatrix.print("PBC_spinMatrix:");
}