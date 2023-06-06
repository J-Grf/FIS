#ifndef _GM_INCLUDE_
#define _GM_INCLUDE_

#include <utility>

#include "matrix.hpp"
#include "preconditioner.hpp"

template<typename T>
using matrixType = std::vector<std::vector<T>>;

std::pair<matrixType<double>, matrixType<double>> gramSchmidt(const Matrix& A, const std::vector<double>& r0, const size_t m);

std::vector<double> getKrylov(const Matrix& A, matrixType<double>& V, matrixType<double>& H, const size_t j, 
const PreConditioner PreCon, std::unique_ptr<ILUout>& ILUobj);

void printKrylov(const matrixType<double>& v);

void printGM(std::pair<matrixType<double>, matrixType<double>> res, const size_t m);

inline bool checkOrthogonality(const std::vector<double> a, const std::vector<double> b) {
    std::cout << "dotP: " << dotP(a, b) << std::endl;
    if (abs(dotP(a, b)) < 10 * eps<double>) 
        return true;
    else 
        return false; 
}

#endif