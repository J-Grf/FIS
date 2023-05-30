#ifndef _PRE_INCLUDE_
#define _PRE_INCLUDE_

#include "matrix.hpp"

enum PreConditioner {
    NONE,
    JACOBI,
    GAUSSSEIDEL,
    ILU
};

void applyPreConditioner(const Matrix& A, std::vector<double>& x, const PreConditioner PreCon);
Matrix Ilu(const Matrix& A);

#endif