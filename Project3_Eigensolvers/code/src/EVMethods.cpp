#include "EVMethods.hpp"

#include <cassert>
#include <unordered_map>
#include <chrono>

PowItObj::PowItObj() {
    diff.resize(bufSize);
    times.resize(bufSize);
    error.resize(bufSize);
}

double powerIteration(const Matrix& A, const std::vector<double>& q0, std::unique_ptr<PowItObj>& PowItPtr, double eps) {
    assert(A.getDim() == q0.size());
    using namespace std::chrono;

#ifndef DISABLEIO
    if(PowItPtr == nullptr) {
        throw std::string("PowItPtr is nullptr, please use DISABLEIO flag when compiling!");
    }
    std::cout << "Using tolerance eps = " << eps << " for power iteration" << std::endl;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
#endif

    std::vector<double> z(A.getDim()), q(q0), Aq(vectorProduct(A, q0));
    double lambdaMax = 0.0, lambdaOld = 0.0, diff = 0.0;
    size_t k = 0;
    while(diff > eps || k < 2) {
        double tmp = norm2(Aq);
        q = Aq / tmp;
        Aq = vectorProduct(A, q);
        lambdaMax = dotP(q, Aq);
        diff = abs(lambdaMax - lambdaOld);

        //To plot error over runtime
#ifndef DISABLEIO
        PowItPtr->diff[k] = diff;
        PowItPtr->error[k] = abs(lambdaCG - lambdaMax);
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        auto time = duration_cast<nanoseconds>(t2 - t1);
        PowItPtr->times[k] = time.count() * 1e-9;
        std::cout << k << " " << PowItPtr->diff[k] << " " << std::setprecision(17) << std::scientific << lambdaMax << std::endl;
#endif
        lambdaOld = lambdaMax;
        k++;
    }

#ifndef DISABLEIO
    PowItPtr->lambdaMax = lambdaMax;
    PowItPtr->maxIdx = k;
#endif

    return lambdaMax;
}

double LanczosMethod(const Matrix& A, std::vector<double>& v, const size_t m, double customEps) {
    
    const std::unordered_map<size_t, double> tolMap = {{3, 1E-1}, {30, 1E-2}, {50, 1E-4}, {75, 1E-6}, {100, 1E-10}};

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
    //length of beta is len(alpha)-1
    beta.erase(beta.cend()-1);

#ifndef DISABLEIO
    printAlphaBeta(alpha, beta);
#endif

    //To MSR format -> m-diagonal elements + 0 to mark start of off-diagonals in VM
    arrayS += m + 1;
    Matrix TM(m, arrayS, 's');
    size_t j = m + 1; 
    for(size_t i = 0; i < TM.getDim(); i++){
        TM.a_VM(i) = alpha[i];
        
        //if off-diagonal element is not 0, add
        if(i < TM.getDim() - 1){
            if(abs(beta[i]) > EPS){
                TM.a_VM(j) = beta[i];
            }
            TM.a_JM(i) = j + 1;
            TM.a_JM(j) = i + 2;
            j++;
        }
    }
    TM.a_JM(TM.getDim() - 1) = arrayS + 1;
    TM.a_JM(TM.getDim()) = arrayS + 1;

    /*std::cout << "-----TM Matrix-------" << std::endl;
    TM.print();*/

    const std::vector<double> q0(m, 1/sqrt(m));
    
    std::unique_ptr<PowItObj> dummy = nullptr;
    
    double _eps = 0.0;
    if(customEps > 0.0) {
        _eps = customEps;
    } else {
        _eps = tolMap.at(m);
    }
    const double res = powerIteration(TM, q0, dummy, _eps);

    return res;
}

void printAlphaBeta(const std::vector<double>& alpha, const std::vector<double>& beta){
    std::string a = "[";
    for(size_t i = 0; i < alpha.size(); i++) {
        a += std::to_string(alpha[i]) + ", ";
    }
    a += " ]";
    std::string b = "[";
    for(size_t i = 0; i < beta.size(); i++) {
        b += std::to_string(beta[i]) + ", ";
    }
    b += " ]";

    std::cout << "alpha: " << a << std::endl    
              << "beta: " << b << std::endl;
}