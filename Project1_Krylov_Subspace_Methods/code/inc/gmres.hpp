#ifndef _GMRES_INCLUDE_
#define _GMRES_INCLUDE_

#include "matrix.hpp"

template<typename T>
std::vector<T> backwardSub(const MatrixCoo& A, const std::vector<T>& b, const size_t m) {

    //use strict upper matrix of Hessenberg matrix
    assert(A.n == b.size() - 1 && "Dimensions of A and b do not coincide, backwardSub not possible!");
    assert(A.n == A.m - 1 && "A is not a square matrix, backwardSub not possible!");
    
    std::vector<T> x(m, 0);
    size_t idx = m - 1;
    std::vector<T> diag = A.getDiagonals();
    
    for(size_t i = 0; i < diag.size(); i++)
        std::cout << "diag[" << i << "]: " << diag[i] << std::endl;


    //first element for backward Sub
    if(diag.size() < A.n) {
        std::cerr << "zeros on diagonal!" << std::endl;
    }
    
    x[idx] = b[idx] / diag[idx];
    for (size_t i = A.values.size() - 3; i --> 0;) {
        const size_t j = A.rows[i];
        const size_t k = A.columns[i];
        
        //ignore elements below diagonal
        if(j > k)
            continue;
        
        x[j] = b[j];
        // substract off-diagonals
        for (size_t n = i; n < A.values.size(); n++) {
            const size_t j_inner = A.rows[n];
            const size_t k_inner = A.columns[n];
            if(k_inner > j_inner)
                x[j_inner] -= A.values[n] * x[k_inner];
        }           
        x[j] /= diag[j];
    } 

    return x; 
}

template<typename T>
std::vector<T> MR_method(const Matrix& A, const std::vector<T>& b, const std::vector<T>& x0) {
    using namespace std;

    vector<T> r0 = b - vectorProduct(A, x0);
    vector<T> p,r,x(A.getDim());
    r = r0;
    T alpha;

    while( norm2(r)/norm2(r0) > Eps) {

        p = vectorProduct(A, r);
        alpha = dotP(r, p) / dotP(p, p);
        x += alpha * r;
        r -= alpha * p;
    } 

    return x;
}

template<typename T>
std::pair<std::vector<T>, T> GMRES(const Matrix& A, const std::vector<T>& x, const std::vector<T>& r0, const size_t m) {
    using namespace std;
    const size_t num = m + 1;

    matrixType<T> v(1);
    const T InvNormR0 = 1.0 / norm2(r0); 
    
    // H: upper Hessenberg 
    MatrixCoo H;

    // v1
    for(size_t i = 0; i < r0.size(); i++) {
        v.at(0).push_back(InvNormR0 * r0.at(i));
    }

    vector<T> e1(num,0); e1[0] = 1;
    vector<T> g = norm2(r0) * e1;

    vector<T> c, s, hj;
    size_t j;
    for(j = 0; j < m; j++) {

        std::cout << "-------GMRES iteration " << j << " --------" <<  std::endl;
        hj = getKrylov(A, v, H, j);
        H.detDimensions();
        H.print();
        for(size_t i = 0; i < hj.size(); i++) {
            std::cout << "hj[" << i << "]: " << hj[i] << std::endl;
        }
        printKrylov(v,m);

        // check mark --- works until here

        // happy break down occured, Krylov space will become linearly dependent
        if(v.size() == j + 1)
            break;

        //Apply previous rotations
        T tmp;
        for(size_t k = 1; k < j; k++) {
            tmp = hj[k - 1]; // store before overwriting
            hj[k - 1] = c[k-1] * tmp + s[k-1] * hj[k];
            hj[k] = -s[k-1] * tmp + c[k-1] * hj[k];
        }

        //new rotations
        const T alpha = sqrt(hj[j] * hj[j] + hj[j+1] * hj[j+1]);
        c.push_back(hj[j] / alpha);
        s.push_back(hj[j+1] / alpha);
        hj[j] = alpha;
        H.setDiagonal(j) = hj[j];

        g[j+1] = -s[j] * g[j];
        g[j] = c[j] * g[j];
        H.print();
    }

    size_t m_tilde = min(j + 1 , m);
    std::cout << "m_tilde: " << m_tilde << std::endl;
    std::vector<T> xm(x.size()), y;

    // back ward substitution
    y = backwardSub(H, g, m_tilde);
    for(size_t i = x.size(); i --> m_tilde; ) {
        xm[i] = x[i];
    }
    std::vector<T> tmp = VP(v,y);
    for(size_t i = 0; i < y.size(); i++) {
        xm[i] = x[i] + tmp[i];
    }

    T rho = abs(g[m_tilde]);
    return make_pair(xm, rho);
}

template<typename T>
std::vector<T> GMRES_Res(const Matrix& A, const std::vector<T>& x0, const std::vector<T>& b, const size_t m){
    assert(A.getDim() == b.size());
    using namespace std;
    vector<T> r0 = b - vectorProduct(A, x0);
    const T r0Norm  = norm2(r0);

    size_t it = 0;
    T rho = 1.0;
    vector<T> x = x0;
    vector<T> r;
    while(rho > Eps) {
        it++;

        std::cout << "-------Iteration number " << it << std::endl;
        const pair<vector<T>, T> res = GMRES(A, x, r0, m);
        x = res.first;
        rho = res.second / r0Norm;

        // save residual
        r.push_back(res.second);

        std::cout << "residual " << r.back() << std::endl;
        std::cout << "rel residual " << rho << std::endl;
        
        for(size_t i = 0; i < x.size(); i++){
            std::cout << "x[" << i << "]: " << x[i] << std::endl;
        }
        
        exit(0);
    }

    return x;
}

#endif
