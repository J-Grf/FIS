#ifndef _GMRES_INCLUDE_
#define _GMRES_INCLUDE_

#include "matrix.hpp"
#include "gram_schmidt.hpp"
#include "preconditioner.hpp"

std::vector<double> backwardSub(const MatrixCoo& A, const std::vector<double>& b, const size_t m);

std::vector<double> backwardSub(const matrixType<double>& A, const std::vector<double>& b, const int m);

std::vector<double> MR_method(const Matrix& A, const std::vector<double>& b, const std::vector<double>& x0);

std::pair<std::vector<double>, double> GMRES(const Matrix& A, const std::vector<double>& x0, const std::vector<double>& b, const size_t m, 
const PreConditioner PreCon);

std::vector<double> GMRES_Res(const Matrix& A, const std::vector<double>& x0, const std::vector<double>& b, const size_t m, const PreConditioner PreCon);

void saveRelResiduals(const std::vector<double>& relRes, const PreConditioner PreCon);
void saveDotPofKrylovVectors(const matrixType<double>& V);

#endif