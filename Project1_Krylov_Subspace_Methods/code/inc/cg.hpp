#ifndef _CG_INCLUDE_
#define _CG_INCLUDE_

#include "matrix.hpp"

std::pair<std::vector<double>, std::vector<double>> CGMethod(const Matrix& A, const std::vector<double>& b, 
const std::vector<double>& x0, const size_t m);

void printCG(const std::vector<double>& x, const std::vector<double>& r, const std::vector<double>& p, double relRes, const size_t it);

void saveData(const std::vector<double>& eANorm, const std::vector<double>& rNorm);

#endif
