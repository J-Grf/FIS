#include "preconditoner.hpp"

Matrix Ilu(const Matrix& A) {
    Matrix P = A;
    
    const size_t n = A.dim;
    const size_t nz = A.array_size;
    std::vector<int> JR(n), JD(n); //PJM = A.JM;
    //std::vector<double> PVM = A.VM;

    //check nochmal mit index shift;
    for(size_t i = 0; i < n; i++){
        JR[i] = i;
        size_t j = 0;
        for(j = A.JM[i]; j < A.JM[i+1] - 1; j++){
            int jc = A.JM[j];
            JR[jc] = j;
            if (jc > i && JD[i] == 0) {
                JD[i] = j;
            }
        }

        if(JD[i] == 0){
            JD[i] = j;
        }
        for(j = A.JM[i]; j < JD[i] - 1; j++) {
            size_t jc = A.JM[j]; 
            P.VM[j] /= P.VM[jc];
            for(size_t jj = JD[jc]; jj < A.JM[jc + 1] - 1; jj++) {
                size_t jk = JR[A.JM[jj]];
                if(jk != 0) {
                    P.VM[jk] -= P.VM[j] * P.VM[jj];
                }
            }
        }

        JR[i] = 0;
        for(j = A.JM[i]; j < A.JM[i+1] - 1; j++) {
            JR[A.JM[j]] = 0;
        }
    }

    return P;
}

/*
    Apply preconditioning:
        - Jacobi: M = D
        - Gauss Seidel: M = L_A (lower triangle of  A)
        - ILU incomplete LU factorization
*/
void applyPreConditioner(const Matrix& A, std::vector<double>& x, const PreConditioner PreCon) {
    switch(PreCon) {
        case JACOBI: {
            for(size_t i = 0; i < A.getDim(); i++) {
                x[i] *= 1.0 / A.getVM(i);
            }
            break;
        }
        case GAUSSSEIDEL: {
            std::vector<double> btmp = x;
            x = forwardSubMSR(A, btmp, A.getDim());
            break;
        }
        case ILU: {
        }
            break;
    } 
}