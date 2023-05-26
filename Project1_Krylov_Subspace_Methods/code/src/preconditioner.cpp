#include "preconditoner.hpp"

Matrix JacobiPre(const Matrix& A) {
    Matrix InvM;

    InvM.dim = A.dim;
    //copy diagonals
    for(size_t i = 0; i < A.dim; i++) {
        InvM.VM.push_back(1.0 / A.VM[i]);
    } 

    return InvM;
}

void applyPreConditioner(const Matrix& A, std::vector<double>& x, const PreConditioner PreCon) {
    Matrix InvM;
    switch(PreCon) {
        case JACOBI:
            InvM = JacobiPre(A);
            for(size_t i = 0; i < A.getDim(); i++) {
                x[i] *= InvM.getVM(i);
            }
            break;
        case GAUSSSEIDEL:
            break;
        case ILU:
            break;
    } 
}