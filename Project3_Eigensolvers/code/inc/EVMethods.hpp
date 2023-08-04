#ifndef _EVMETHODS_INCLUDE_
#define _EVMETHODS_INCLUDE_

#include <vector>
#include "matrix.hpp"

const double EPS = 1E-12;

double powerIteration(const Matrix& A, const std::vector<double>& q0, double eps = 1E-8);

double LanczosMethod(const Matrix& A, std::vector<double>& v, const size_t m);

#endif