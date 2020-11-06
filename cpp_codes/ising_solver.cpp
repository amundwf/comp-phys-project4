// The solver class for running the Ising model.

#include "ising_solver.hpp"

using namespace std;
using namespace arma;

// Constructor with an initial spin matrix as input:
IsingSolver::IsingSolver(imat spinMatrix){
    this->spinMatrix = spinMatrix;
}

void IsingSolver::printSpinSystem(){
    // This function prints the current spin system (the spinSystem
    // matrix).
    spinMatrix.print("spinMatrix:");
}