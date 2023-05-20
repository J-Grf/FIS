#ifndef _GMRES_INCLUDE_
#define _GMRES_INCLUDE_

#include "matrix.hpp"

template<typename T>
std::vector<T> MR_method(const Matrix& A, const std::vector<T>& b, const std::vector<T>& x0) {
    using namespace std;

    vector<T> r0 = b - vectorProduct(A, x0);
    vector<T> p,r,x(A.getDim());
    r = r0;
    T alpha = 0.0;

    do {
        p = vectorProduct(A, r);
        alpha = dotP(r, p) / dotP(p, p);
        x += alpha * r;
        r -= alpha * p;
    } while (norm2(r) < eps<T>);

    return x;
}

#endif
