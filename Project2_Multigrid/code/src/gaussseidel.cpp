#include "gaussseidel.hpp"
#include <iostream>

m_type GaussSeidel(m_type u, m_type f, const size_t nu) {
    const size_t N = u.size();
    const double hsq = pow(1 / N ,2);

    /*
        - u[i-1][j] and u[i][j-1] at k
        - u[i+1][j] and u[i][j+1] at k - 1 from previous iteration!
    */
    m_type u_old = u;
    double infNorm = -1, diffU = -1;
    for(size_t k = 0; k < nu; k++) {
        for(size_t i = 1; i < N - 1; i++) {
            for(size_t j = 1; j < N - 1; j++) {
                u[i][j] = 0.25 * (hsq * f[i][j] + u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1]);
                diffU = abs(u[i][j] - u_old[i][j]);
                if(diffU > infNorm) {
                    infNorm = diffU;
                } 
            }
        }
        u_old = u;
        std::cout << " ||u_nu - u_{nu-1}|| " << infNorm << std::endl;
    }

    return u;
}