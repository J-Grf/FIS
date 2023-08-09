#ifndef _EVMETHODS_INCLUDE_
#define _EVMETHODS_INCLUDE_

#include <vector>
#include "matrix.hpp"

constexpr double EPS = 1E-12;
constexpr double lambdaCG = 9.5986080894852857E3;

double powerIteration(const Matrix& A, const std::vector<double>& q0, double eps = 1E-8);

double LanczosMethod(const Matrix& A, std::vector<double>& v, const size_t m);

void printAlphaBeta(const std::vector<double>& alpha, const std::vector<double>& beta);

#endif