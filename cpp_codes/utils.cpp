#include "utils.hpp"
#include "ising_solver.hpp"
#include <iostream>
#include <cmath>
#include <armadillo>
#include <fstream>

using namespace std;
using namespace arma;

void writeIntegerMatrixToCSV(imat spinMatrix, string filename, string directory){
    // This function saves an integer matrix (imat), e.g. the spin matrix, in a .csv file.
    ofstream ofile;
    string filePath = directory + filename;
    // Save the matrix in CSV format in the specified directory:
    spinMatrix.save(csv_name(filePath));
}


void run_4c_ising(){
    int L = 2; // 2x2 spin system
    double J = 1;
    IsingSolver isingSolver2x2(L, J); // Generates a random 2x2 spin state.

}