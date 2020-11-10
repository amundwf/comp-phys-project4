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

void unit_testing_2x2(){
    // This function reproduces the table for the 2x2 spin system (energy,
    // magnetization and number of up spins).
    // Test writing some binary (values -1 and 1) spin matrix (2x2) to a .csv file:
    int L = 2; double J = 1;
    imat spin_mat(L, L, fill::ones);

    // For all 16 different possible spin configurations:
    int num = 0;
    int n_up = 0;
    for (int s1=-1; s1<=1; s1+=2){
        for (int s2=-1; s2<=1; s2+=2){
            for (int s3=-1; s3<=1; s3+=2){
                for (int s4=-1; s4<=1; s4+=2){
                    n_up = 0;
                    spin_mat(0,0) = s1; spin_mat(0,1) = s2;
                    spin_mat(1,0) = s3; spin_mat(1,1) = s4;
                    if (s1==1) n_up += 1; if (s2==1) n_up += 1;
                    if (s3==1) n_up += 1; if (s4==1) n_up += 1;
                    

                    // Make an ising solver object:
                    IsingSolver isingSolver2x2(spin_mat, J);

                    // Calculate energy and magnetization:
                    double E = isingSolver2x2.calculate_E();
                    int M = isingSolver2x2.calculate_M();

                    num += 1;
                    cout << "spin_mat " << to_string(num) << endl;
                    // Uncomment this if you want to show the spin matrix:
                    //spin_mat.print();

                    cout << "number of up spins: " << n_up;
                    if (n_up==2){if (s1==s4) cout << " (diagonal)"; else cout << " (non-diagonal)";}

                    cout << "\nEnergy: " << E << endl;
                    cout << "Magnetization: " << M << "\n\n";
                }
            }
        }
    }
}

void run_4c_ising(){
    int L = 2; // 2x2 spin system
    double J = 1;
    IsingSolver isingSolver2x2(L, J); // Generates a random 2x2 spin state.

}