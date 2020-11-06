#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <cmath>
#include <armadillo>
#include <fstream>

void writeIntegerMatrixToCSV(arma::imat spinMatrix, std::string filename, std::string directory);

void run_4c_ising();


#endif