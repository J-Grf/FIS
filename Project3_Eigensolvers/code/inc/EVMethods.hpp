#ifndef _EVMETHODS_INCLUDE_
#define _EVMETHODS_INCLUDE_

#include <vector>
#include <cassert>
#include <map>
#include "matrix.hpp"

const double EPS = 1E-12;

double powerIteration(const Matrix& A, const std::vector<double>& q0, double eps = 1E-8) {
    assert(A.getDim() == q0.size());
    
    std::vector<double> z(A.getDim()), q(q0);
    double lambdaMax = 0.0, lambdaOld = 0.0, diff = 0.0;
    while(diff < eps) {
        z = vectorProduct(A, q);
        double tmp = norm2(z);
        q = z / tmp;
        lambdaMax = dotP(q, vectorProduct(A, q));
        diff = abs(lambdaMax - lambdaOld);
        lambdaOld = lambdaMax;
    }

    return lambdaMax;
}

double LanczosMethod(const Matrix& A, std::vector<double>& v, const size_t m) {
    
    const std::map<size_t, double> tolMap = {{30, 1E-2}, {50, 1E-4}, {75, 1E-6}, {100, 1E-10}};

    std::vector<double> w, beta(m), alpha(m), vOld(v);
    size_t arrayS = 0;
    for(size_t j = 0; j < m; j++){
        if(j == 0) {
            w = vectorProduct(A, v);
        } else {
            w = vectorProduct(A, v) - beta[j-1] * vOld;
        }
        alpha[j] = dotP(v, w);
        w -= alpha[j] * v;
        beta[j] = norm2(w);
        vOld = v;
        v = w / beta[j];
        if(abs(beta[j]) > EPS) 
            arrayS++;
    }

    //To MSR format -> m-diagonal elements + 0 to mark start of off-diagonals in VM
    arrayS += m + 1;
    Matrix TM(m, arrayS, 's');
    size_t j = m + 1; 
    //TODO: Test
    for(size_t i = 0; i < TM.getDim(); i++){
        TM.a_VM(i) = alpha[i];
        
        //if off-diagonal element is not 0, add
        if(abs(beta[i]) > EPS)
            TM.a_VM(j) = beta[i];
            TM.a_JM(i) = j + 1;
            TM.a_JM(j) = i;
            j++;
    }

    const std::vector<double> q0(1/sqrt(m), m);
    
    return powerIteration(TM, q0, tolMap[m]);
}


#endif