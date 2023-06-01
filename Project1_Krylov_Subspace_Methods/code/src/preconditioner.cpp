#include "preconditioner.hpp"
#include "matrix.hpp"

ILUout::ILUout(const size_t _nz, const size_t _n) : nz(_nz), n(_n), PV(std::vector<double>(_nz)), 
PJM(std::vector<int>(_nz)), JD(std::vector<int>(_n)) {}

ILUout Ilu(const Matrix& A) {
    const size_t n = A.getDim();
    const size_t nz = A.getArrSize();

    ILUout P(nz, n);
    P.PJM = A.getJMVec();
    P.PV = A.getVMVec();
    std::vector<int> JR(n), JD(n);

    //check nochmal mit index shift;
    for(size_t i = 0; i < n; i++){
        JR[i] = i;
        size_t j = 0;
        for(j = A.a_JM(i); j < A.a_JM(i+1) - 1; j++){
            int jc = A.a_JM(j);
            JR[jc] = j;
            if (jc > i && JD[i] == 0) {
                JD[i] = j;
            }
        }

        if(JD[i] == 0){
            JD[i] = j;
        }
        for(j = A.a_JM(i); j < JD[i] - 1; j++) {
            size_t jc = A.a_JM(j); 
            P.PV[j] /= P.PV[jc];
            for(size_t jj = JD[jc]; jj < A.a_JM(jc + 1) - 1; jj++) {
                size_t jk = JR[A.a_JM(jj)];
                if(jk != 0) {
                    P.PV[jk] -= P.PV[j] * P.PV[jj];
                }
            }
        }

        JR[i] = 0;
        for(j = A.a_JM(i); j < A.a_JM(i+1) - 1; j++) {
            JR[A.a_JM(j)] = 0;
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
                x[i] *= 1.0 / A.a_VM(i);
            }
            break;
        }
        case GAUSSSEIDEL: {
            std::vector<double> btmp = x;
            x = forwardSubMSR(A, btmp, A.getDim());
            break;
        }
        case ILU: {
            ILUout M = Ilu(A);
            Matrix M_dummy = A;
            M_dummy.getJMVec() = M.PJM;
            M_dummy.getVMVec() = M.PV;

            std::vector<double> btmp = x;
            std::vector<double> y = forwardSubMSR(M_dummy, btmp, M_dummy.getDim());
            x = backwardSubMSR(M_dummy, y, M_dummy.getDim());
            break;
        }
        default: {
            std::cout << "preconditioner unknown!" << std::endl;
            break;
        }
    } 
}