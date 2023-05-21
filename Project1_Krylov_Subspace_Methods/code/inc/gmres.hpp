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

template<typename T>
std::pair<std::vector<T>, T> GMRES(const Matrix& A, const std::vector<T>& x, const std::vector<T>& r0, const T m) {
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

    //TODO, check if m is correct here
    vector<T> e1(m,0); e1[0] = 1;
    vector<T> g = norm2(r0) * e1;

    vector<T> c, s;
    for(size_t j = 0; j < m; j++) {
        getKrylov(A, v, H, j);
        for(size_t k = 1; k < j; k++) {
            H[k - 1][j] = c[k-1] * H[k-1][j] + s[k-1] * H[k][j];
            H[k][j] = -s[k-1] * H[k-1][j] + c[k-1] * H[k][j];
        }
        const T tmp = sqrt(H[j][j] * H[j][j] + H[j+1][j] * H[j+1][j]);
        c.push_back(H[j][j] / tmp);
        s.push_back(H[j+1][j] / tmp);
        H[j][j] = c[j] * H[j][j] + s[j] * H[j+1][j];

        g[j+1] = -s[j] * g[j]; //add exit condition
        g[j] = c[j] * g[j];
    }

    //TODO compute x
    std::vector<T> xm;
    T rho = abs(g[m]);

    return make_pair(xm, rho);
}

template<typename T>
std::vector<T> GMRES_Res(const Matrix& A, const std::vector<T>& x0, const std::vector<T>& b, const T m){
    assert(A.getDim() == b.size());
    using namespace std;
    vector<T> r = b - vectorProduct(A, x0);

    T rho  = norm2(r);
    vector<T> x = x0;
    while(rho > eps<T>) {
        const pair<vector<T>, T> res = GMRES(A, x, r, m);
        x = res.first;
        rho = res.second;
    }

    return x;
}

#endif
