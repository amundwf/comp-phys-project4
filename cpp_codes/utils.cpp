#include "utils.hpp"
#include "ising_solver.hpp"
#include "omp.h"
#include <time.h>
#include <iostream>
#include <cmath>
#include <armadillo>
#include <fstream>
#include <random>
#include <iomanip>
#define NUM_THREADS 8

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
    int N_MC = 1e4;//
    // Choose the temperature value to run the Metropolis algorithm for:
    double T = 1;

    imat spin_mat(L, L, fill::ones);
    spin_mat(0,0) = 1; spin_mat(0,1) = 1;
    spin_mat(1,0) = 1; spin_mat(1,1) = 1;

    // Make the Ising solver object (random 2x2 spin matrix):
    //IsingSolver isingSolver2x2(L, T, N_MC); // Generates a random 2x2 spin state.
    IsingSolver isingSolver2x2(spin_mat, T, N_MC);

    cout << "State before Metropolis:\n";
    isingSolver2x2.print_spinMatrix();

    cout << "\nPerforming Metropolis algorithm...\n\n";
    isingSolver2x2.run_metropolis_full(); // Run Metropolis
    
    cout << "State after Metropolis:\n";
    isingSolver2x2.print_spinMatrix();

    // Results matrix containing E_list, M_list, E2_list, M2_list, M_abs_list:
    mat results_E_M = isingSolver2x2.get_E_list_M_list();

    // Mean values results matrix, containing <E>, <M>, C_V, chi:
    mat results_MVs = isingSolver2x2.get_mean_results();

    // Save the results:
    // The results directory (the saving location of the results):
    string directory = "../results/4c_ising/";
    // First, save the results for E and M (and E^2, M^2 and |M|) as a function
    // of MC cycles:
    string filename = "2x2_T=1.0_N=" + to_string(N_MC) +".csv";
    field<string> header(results_E_M.n_cols);
    header(0) = "MC_cycle"; header(1) = "E"; header(2) = "E^2";
    header(3) = "M"; header(4) = "M_abs"; header(5) = "M^2";
    writeGeneralMatrixToCSV(results_E_M, header, filename, directory);

    // Now save the results of the mean value quantities:
    string filenameMVs = "4c_MVs.csv";
    field<string> headerMVs(results_MVs.n_cols);
    headerMVs(0) = "E_mean"; headerMVs(1) = "E2_mean";
    headerMVs(2) = "M_mean"; headerMVs(3) = "M_abs_mean"; headerMVs(4) = "M2_mean";
    headerMVs(5) = "C_V"; headerMVs(6) = "chi";
    writeGeneralMatrixToCSV(results_MVs, headerMVs, filenameMVs, directory);
}

void run_4d_most_likely_state(){
    int L = 20; // 2x2 spin system
    // (Just some random test values at first when testing):
    int N_MC = 100000;
    // Choose the temperature value to run the Metropolis algorithm for:
    double T = 2.4;

    // Make the Ising solver object (random 20x20 spin matrix):
    //IsingSolver isingSolver20x20(L, T, N_MC); // Generates a random LxL spin state.
    imat spin_mat(L, L, fill::ones);
    IsingSolver isingSolver20x20(spin_mat, T, N_MC);

    // Run the metropolis algorithm, and for each MC cycle, calculate mean values
    // and add them to a list. Then save the lists for plotting.
    // Basically do the same thing as run_metropolis_full(), but now with updating
    // and saving the mean values after each MC cycle (takes longer time).
    // Create array to contain the mean values for all MC cycles:
    mat mean_value_results = zeros(N_MC,8);
    // ^ Columns: MC_cycles, E_mean, E2_mean, M_mean, M_abs_mean, M2_mean, C_V, chi.
    cout << "\nPerforming Metropolis algorithm...\n\n";
    for (int n=1; n<=N_MC; n++){ // N_MC MC cycles
        // Run one MC cycle:
        isingSolver20x20.metropolis_one_time_and_update_E_M_lists(n);

        // Calculate and update the mean values:
        isingSolver20x20.calculate_mean_results_current_MC_cycle(n);

        // Update the mean results array:
        mean_value_results(n-1,0) = n; // MC cycle number.
        mean_value_results(n-1, span(1,7)) = isingSolver20x20.get_mean_results();
    }   
    //cout << "State after Metropolis:\n";
    //isingSolver20x20.print_spinMatrix();

    // Save the results:
    // The results directory (the saving location of the results):
    string directory = "../results/4d/";
    // First, save the results for E and M (and E^2, M^2 and |M|) as a function
    // of MC cycles:
    string filename = "4d_mean_values_vs_" + to_string(N_MC) + "_cycles_T=2.4_ordered.csv";
    field<string> header_MVs(mean_value_results.n_cols);
    header_MVs(0) = "MC_cycle"; header_MVs(1) = "E_mean"; header_MVs(2) = "E2_mean";
    header_MVs(3) = "M_mean"; header_MVs(4) = "M_abs_mean"; header_MVs(5) = "M2_mean";
    header_MVs(6) = "C_V"; header_MVs(7) = "chi";
    writeGeneralMatrixToCSV(mean_value_results, header_MVs, filename, directory);

    // Save E_list, E2_list, etc. in another file:
    // Results matrix containing E_list, M_list, E2_list, M2_list, M_abs_list:
    mat results_E_M = isingSolver20x20.get_E_list_M_list();
    // The results directory (the saving location of the results):
    // First, save the results for E and M (and E^2, M^2 and |M|) as a function
    // of MC cycles:
    filename = "4d_results_"+ to_string(N_MC) +"_T=2.4_ordered.csv";
    field<string> header(results_E_M.n_cols);
    header(0) = "MC_cycle"; header(1) = "E"; header(2) = "E^2";
    header(3) = "M"; header(4) = "M_abs"; header(5) = "M^2"; header(6) = "flipsAccepted";
    writeGeneralMatrixToCSV(results_E_M, header, filename, directory);
}


void run_4f(){

    // Monte Carlo cycles.
    int N_MC = 100000;

    // Set the number of threads.
    omp_set_num_threads(NUM_THREADS);
    cout << "The number of processors available = " << omp_get_num_procs ( ) << endl;
    
    // Set the temperature steps.
    double Tstart = 2.1;
    double Tend = 2.4;
    double Tsteps = 32;
    double T_step_length = (Tend - Tstart)/Tsteps;

    // Make data into strings.
        // Create an output string stream
    std::ostringstream streamObj3;
    // Set Fixed -Point Notation
    streamObj3 << std::fixed;
    // Set precision to 1 digits
    streamObj3 << std::setprecision(1);
    //Add double to stream
    streamObj3 << Tstart;
    std::string stringTstart = streamObj3.str();

    std::ostringstream streamObj;
    // Set Fixed -Point Notation
    streamObj << std::fixed;
    // Set precision to 1 digits
    streamObj << std::setprecision(1);
    streamObj << Tend;
    std::string stringTend = streamObj.str();

    // Temperature vector.
    vector<double> T_schedule;
    for (double i = Tstart; i <= Tend; i+= T_step_length){
        T_schedule.push_back(i);
    }

    // Start timing.
    double wtime = omp_get_wtime ( );

    // Needs to be declared before parallel task.
    int i;
    double T;

    // Loop over spin matrices.
    for (int L = 40; L <= 100; L+=20){
        // Matrix for results.
        mat mean_value_results = zeros(Tsteps,8);
        cout << "The spin matrix dimentions is: " << L << "x" << L << endl;
        cout << "Running Metropolis algo ..." << endl;
        // Start parallel task.
        # pragma omp parallel shared(mean_value_results, L, N_MC, Tsteps, T_schedule) private (i)
        {
        # pragma omp for
        for (int i=0; i < int(Tsteps); i++){

            double T = T_schedule[i];
            //cout << "\nT is: " << T << endl;
            IsingSolver isingSolver(L, T, N_MC);
            
            // Begin equilibrium cycles
            for (int n=1; n<=int(N_MC*0.1); n++){ 
                isingSolver.metropolis_one_time_parallel();
            }  
            // Clear data.
            isingSolver.init_parallel_variables();
            
            // Begin recording real data.
            for (int n=int(N_MC*0.1); n<=N_MC; n++){ 
                isingSolver.metropolis_one_time_parallel();
            }  
            mean_value_results(i, 0) = T;
            mean_value_results(i, span(1,7)) = isingSolver.get_mean_results_parallel();
        } 
        } // End parallel region. only one processor for I/O.
        // Print wall time.
        wtime = omp_get_wtime ( ) - wtime;
        cout << "  Elapsed time in seconds = " << wtime << endl;
        
        // Save the results.
        string directory = "../results/4f/";
        string filename = stringTstart + "_" + stringTend + "_steps=" + to_string(int(Tsteps)) + "_L=" + to_string(L) + "_N=" + to_string(N_MC) + ".csv";
        field<string> header(mean_value_results.n_cols);
        header(0) = "Temperature"; header(1) = "E_mean"; header(2) = "E2_mean"; header(3) = "M_mean";
        header(4) = "M_abs_mean"; header(5) = "M2_mean"; header(6) = "C_V";
        header(7) = "chi"; 
        writeGeneralMatrixToCSV(mean_value_results, header, filename, directory);
    }
//end
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