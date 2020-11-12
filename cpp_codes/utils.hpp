#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <cmath>
#include <armadillo>
#include <fstream>

void writeIntegerMatrixToCSV(arma::imat spinMatrix, std::string filename, std::string directory);

int randrange_int(int a, int b);

double randrange_float(double a, double b);

int random_plusMinus1(); // Returns randomly either +1 or -1 (uniformly).

void unit_testing_2x2();

void run_4c_ising();


#endif