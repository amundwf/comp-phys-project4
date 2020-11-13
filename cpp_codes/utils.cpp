#include "utils.hpp"
#include "ising_solver.hpp"
#include <iostream>
#include <cmath>
#include <armadillo>
#include <fstream>
#include <random>

using namespace std;
using namespace arma;

void writeIntegerMatrixToCSV(imat spinMatrix, string filename, string directory){
    // This function saves an integer matrix (imat), e.g. the spin matrix, in a .csv file.
    ofstream ofile;
    string filePath = directory + filename;
    // Save the matrix in CSV format in the specified directory:
    spinMatrix.save(csv_name(filePath));
}

int randrange_int(int a, int b){
    // Returns a random integer between integers a and b, including both endpoints,
    // meaning that this function can also return a and b, not just the numbers in
    // between.

    // Seed the random number generator and get a random number:
    std::random_device rand_dev; // This works kind of like a seed?
    std::mt19937 generator(rand_dev()); // Mersenne twister, a pseudo-random number
    // generator.
    std::uniform_int_distribution<int> distr(a, b);
    // Get a random number:
    int randInt = distr(generator);
    // This gives a 50/50 chance between a positive and a negative number.
    // Therefore, simply return the sign of r (either +1 or -1).
    return randInt;
}

double randrange_float(double a, double b){
    // Returns a random floating point number between numbers a and b.

    // Seed the random number generator and get a random number:
    std::random_device rand_dev; // This works kind of like a seed?
    std::mt19937 generator(rand_dev()); // Mersenne twister, a pseudo-random number
    // generator.
    std::uniform_real_distribution<double> distr(a, b);
    // Get a random number:
    double randNum = distr(generator);
    // This gives a 50/50 chance between a positive and a negative number.
    // Therefore, simply return the sign of r (either +1 or -1).
    return randNum;
}

int random_plusMinus1(){
    double randNum = randrange_float(-1,1);
    // This gives a 50/50 chance between a positive and a negative number.
    // Therefore, simply return the sign of r (either +1 or -1).
    return sign(randNum);
}

void unit_testing_2x2(){
    // This function reproduces the table for the 2x2 spin system (energy,
    // magnetization and number of up spins).
    // Test writing some binary (values -1 and 1) spin matrix (2x2) to a .csv file:
    int L = 2;
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
                    IsingSolver isingSolver2x2(spin_mat, 1, 1);

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
    // (Just some random test values at first when testing):
    int N_MC = 10000; //int N_MC = 10;
    // Choose the temperature value to run the Metropolis algorithm for:
    double T = 1;

    // File name and directory (saving location of results):
    // (deal with this after you get the algorithm itself to work)
    // Save <E>, <M>, C_V, Chi etc. in this file, as function of either
    // temperature or # MC cycles performed? Or is one of the mean values
    // calculated from performing a lot of MC cycles and taking the mean of that,
    // and then that is one mean value for one temperature value? So the results
    // will be functions of temperature? (<-- this makes sense.)
    string filename = "blabla_4c.csv";
    string directory = "../results/4c_ising/";

    // Make the Ising solver object (random 2x2 spin matrix):
    IsingSolver isingSolver2x2(L, T, N_MC); // Generates a random 2x2 spin state.

    // Do the Monte Carlo sampling for the chosen temperature:
    //vec ExpectationValues = zeros<mat>(5);
    //MetropolisSampling(NSpins, MCcycles, Temperature, ExpectationValues);
    // Run the full metropolis algorithm (for N_MC MC cycles):
    //isingSolver2x2.run_metropolis_full();

    cout << "State before Metropolis:\n";
    isingSolver2x2.print_spinMatrix();

    cout << "\nPerforming Metropolis algorithm...\n\n";
    isingSolver2x2.run_metropolis_full(); // Run Metropolis
    
    cout << "State after Metropolis:\n";
    isingSolver2x2.print_spinMatrix();
    //isingSolver2x2.print_E_list_and_M_list();
    //imat results_EVs = isingSolver2x2.get_spinMatrix();

    


    // Results matrix: First column: TList (all values of temperature). Second, third
    // and all subsequent columns: <E>, <M>, C_V, chi, ... (for the corresponding 
    // values of T).
    mat results = isingSolver2x2.get_E_list_and_M_list();

    // Save results matrix, with a header:
    //(Do something like this, from utils.cpp from project 3?: )
    //field<string> header(n_columns);
    //header(0) = "T"; header(1) = "E_mean"; header(2) = "M_mean"; header(3) = "C_V"; header(4) = "chi";
    //writeGeneralMatrixToCSV(results, header, filename, directory);

    field<string> header(results.n_cols);
    header(0) = "MC_cycle"; header(1) = "E_mean"; header(2) = "M_mean"; //header(3) = "C_V"; header(4) = "chi";
    writeGeneralMatrixToCSV(results, header, filename, directory);

}

/* // 4f: Multiple values of T.
void run_4f(){
    int L = 2; // 2x2 spin system
    double J = 1;

    // (Just some random test values at first when testing):
    int N_MC = 10;
    double T_initial = 0.3; double T_final = 1; // Initial and final temperatures.
    double dT = 0.1; // Temperature interval.
    int N_T = round((T_final-T_initial)/dT); // Number of temperature steps.
    // Make a list with the temperature steps:
    vec TList = vec(N_T);
    for (int i=0; i<=N_T-1; i++){TList(i) = T_initial + i*dt;}

    // File name and directory (saving location of results):
    string filename = "blabla_4f.csv";
    string directory = "../results/4f_ising/";

    // Start Monce Carlo sampling by performing a loop over the chosen temperatures:
    for (double i=0; T<=){
        T = TList(i);

        vec ExpectationValues = zeros<mat>(5);
    }

    IsingSolver isingSolver2x2(L, T, N_MC); // Generates a random 2x2 spin state.

    // Results matrix: First column: TList (all values of temperature). Second, third
    // and all subsequent columns: <E>, <M>, C_V, chi, ... (for the corresponding 
    // values of T).
    //mat results = 

    // Save results matrix, with a header:
    //(Do something like this, from utils.cpp from project 3?: )
    //field<string> header(n_columns);
    //header(0) = "T"; header(1) = "E_mean"; header(2) = "M_mean"; header(3) = "C_V"; header(4) = "chi";
    //writeGeneralMatrixToCSV(results, header, filename, directory);

}
*/

void writeGeneralMatrixToCSV(mat results, field<string> columnLabels, string filename, string directory){
    // columnLabels contains the labels of the columns, e.g. "t", "x", "y", "z" or "t", "L".
    // It should have the same amount of elements as 'results' has columns.
    ofstream ofile;
    string filePath = directory + filename;

    // Save matrix in CSV format with the column labels in the header:
    //results.save(csv_name("results.csv", header));
    results.save(csv_name(filePath, columnLabels));
}