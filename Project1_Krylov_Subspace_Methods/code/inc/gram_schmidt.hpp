#ifndef _GM_INCLUDE_
#define _GM_INCLUDE_

#include <utility>
#include <cmath>
#include <cassert>

#include "matrix.hpp"

template<typename T>
using matrixType = std::vector<std::vector<T>>;

/* Gram-Schmidt orthogonalization procedure
   v contains here the different orthonormal vectors per row
   not per column, to allow a faster access in the architecture
*/
template<typename T>
std::pair<matrixType<T>, matrixType<T>> gramSchmidt(const Matrix& A, 
const std::vector<T>& r0, const size_t m) {
    assert(A.getDim() == r0.size());

    using namespace std;

    const size_t num = m + 1;

    matrixType<T> v(1), H;
    const T InvNormR0 = 1.0 / norm2(r0); 
    

    // construct v1
    for(size_t i = 0; i < num; i++) {
        H.push_back(vector<T>(m, 0));
    }
    
    for(size_t i = 0; i < r0.size(); i++) {
        v.at(0).push_back(InvNormR0 * r0.at(i));
    }
    
    std::vector<T> w;
    for(size_t j = 0; j < m; j++) {
        std::cout << "----j: " << j << std::endl;
        w = vectorProduct(A, v[j]);
        for(size_t i = 0; i < j + 1; i++) {
            H.at(i).at(j) = dotP(v.at(i), w);
            w -=  H.at(i).at(j) * v.at(i); 
        }
        H.at(j + 1).at(j) = norm2(w);
        if (H.at(j + 1).at(j) > eps<T>) {
            v.push_back((1.0 / H.at(j + 1).at(j)) * w);
        } else
            break;
    }

    return make_pair(v,H);
}

template<typename T>
std::vector<T> getKrylov(const Matrix& A, matrixType<T>& v, MatrixCoo& H, const size_t j) {
    assert(A.getDim() == v[j].size());

    std::vector<T> w;
    std::vector<T> hj;
    
    std::cout << "----j: " << j << std::endl;
    w = vectorProduct(A, v[j]);

    for(size_t i = 0; i < j + 1; i++) {
        T value = dotP(v.at(i), w);

        //store all non-zero values in sparse format
        if(abs(value) > eps<T>) {
            H.append(i, j, value);
        }

        hj.push_back(value);
        
        w -=  hj.at(i) * v.at(i); 
    }

    T value2 = norm2(w);
    if(value2 > eps<T>) {
        hj.push_back(value2);
        H.append(j + 1, j, value2);
    }
    
    T invNormW = 1.0 / value2;
    if (value2 > eps<T>) {
        v.push_back(invNormW * w);
    } 
    return hj;
}

template<typename T>
void printKrylov(const matrixType<T>& v, const size_t m) {
    
    for(size_t i = 0; i < v.size(); i++) {
        for(size_t j = 0; j < v[0].size(); j++) {
            std::cout << "V[" << i << "][" << j << "]: " << v[i][j] << std::endl;
        }
    }
}

template<typename T>
void printGM(std::pair<matrixType<T>, matrixType<T>> res, const size_t m) {
    std::cout << "res.first[0].size(): " << res.first[0].size() << std::endl
              << "res.second.size(): " << res.second.size() << std::endl;
    //assert(res.first[0].size() == res.second.size());
    for(size_t i = 0; i < m + 1; i++) {
        for(size_t j = 0; j < m; j++) {
            std::cout << "H[" << i << "][" << j << "]: " << res.second.at(i).at(j) << std::endl;
        }
    }
    for(size_t i = 0; i < m + 1; i++) {
        for(size_t j = 0; j < 3; j++) {
            std::cout << "V[" << i << "][" << j << "]: " << res.first.at(i).at(j) << std::endl;
        }
    }
}

template<typename T>
inline bool checkOrthogonality(const std::vector<T> a, const std::vector<T> b) {
    std::cout << "dotP: " << dotP(a, b) << std::endl;
    if (abs(dotP(a, b)) < 10 * std::numeric_limits<T>::epsilon()) 
        return true;
    else 
        return false; 
}

#endif