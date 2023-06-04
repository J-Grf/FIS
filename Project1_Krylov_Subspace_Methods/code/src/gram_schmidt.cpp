#include "gram_schmidt.hpp"

/* Gram-Schmidt orthogonalization procedure
   v contains here the different orthonormal vectors per row
   not per column, to allow a faster access in the architecture
*/
std::pair<matrixType<double>, matrixType<double>> gramSchmidt(const Matrix& A, const std::vector<double>& r0, const size_t m) {
    assert(A.getDim() == r0.size());

    using namespace std;

    const size_t num = m + 1;

    matrixType<double> v(1), H;
    const double InvNormR0 = 1.0 / norm2(r0); 
    
    // construct v1
    for(size_t i = 0; i < num; i++) {
        H.push_back(vector<double>(m, 0));
    }
    
    for(size_t i = 0; i < r0.size(); i++) {
        v[0].push_back(InvNormR0 * r0[i]);
    }
    
    std::vector<double> w;
    for(size_t j = 0; j < m; j++) {
        std::cout << "----j: " << j << std::endl;
        w = vectorProduct(A, v[j]);
        for(size_t i = 0; i < j + 1; i++) {
            H[i][j] = dotP(v[i], w);
            w -=  H[i][j] * v[i]; 
        }
        H[j + 1][j] = norm2(w);
        if (H[j + 1][j] > eps<double>) {
            v.push_back((1.0 / H[j + 1][j]) * w);
        } else
            break;
    }

    return make_pair(v,H);
}

std::vector<double> getKrylov(const Matrix& A, matrixType<double>& V, matrixType<double>& H, const size_t j, const PreConditioner PreCon) {
    assert(A.getDim() == V[j].size());

    std::vector<double> w;
    std::vector<double> hj;

    w = vectorProduct(A, V[j]);
#ifndef DISABLEIO
    if(PreCon != NONE) {
        // w will be modified wbar
        applyPreConditioner(A, w, PreCon);
    }
#endif

    for(size_t i = 0; i < j + 1; i++) {
        double value = dotP(V[i], w);

        hj.push_back(value);
        H[i][j] = value;

        w -=  hj[i] * V[i]; 
    }

    double value2 = norm2(w);
    hj.push_back(value2);
    H[j + 1][j] = value2;
    
    double invNormW = 1.0 / value2;
    V.push_back(invNormW * w);
    
    return hj;
}

void printKrylov(const matrixType<double>& v) {
    for(size_t i = 0; i < v.size(); i++) {
        for(size_t j = 0; j < v[0].size(); j++) {
            std::cout << "V[" << i << "][" << j << "]: " << v[i][j] << std::endl;
        }
    }
}

void printGM(std::pair<matrixType<double>, matrixType<double>> res, const size_t m) {
    for(size_t i = 0; i < m + 1; i++) {
        for(size_t j = 0; j < m; j++) {
            std::cout << "H[" << i << "][" << j << "]: " << res.second[i][j] << std::endl;
        }
    }
    for(size_t i = 0; i < m + 1; i++) {
        for(size_t j = 0; j < 3; j++) {
            std::cout << "V[" << i << "][" << j << "]: " << res.first[i][j] << std::endl;
        }
    }
}