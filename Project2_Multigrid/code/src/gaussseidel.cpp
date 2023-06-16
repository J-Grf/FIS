#include "gaussseidel.hpp"
#include <iostream>

m_type GaussSeidel(m_type& u, const m_type& u_ex, const m_type& f, const size_t nu) {
    const size_t N = u.size();
    const double hsq = pow(1 / static_cast<double>(N) ,2);

    /*
        - u[i-1][j] and u[i][j-1] at k
        - u[i+1][j] and u[i][j+1] at k - 1 from previous iteration!
    */
    m_type u_old = u;
    size_t k = 0;
    double infNorm = 1E10;
    while(infNorm > eps) {
        k++;
        double diffU;
        double max_error = -1, diffUex;
        for(size_t i = 1; i < N - 1; i++) {
            for(size_t j = 1; j < N - 1; j++) {
                u[i][j] = 0.25 * (hsq * f[i][j] + u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1]);
                
                diffU = abs(u[i][j] - u_old[i][j]);
                diffUex = abs(u_ex[i][j] - u[i][j]);
                if(diffU < infNorm) 
                    infNorm = diffU;
                if(diffUex > max_error)
                    max_error = diffUex;
                if (infNorm < eps)
                    break;
            }
        }
        u_old = u;
        std::cout << " ||u_nu - u_{nu-1}|| " << infNorm << " after " << k << " iterations" << std::endl;
        std::cout << "maximum converged error: " << max_error << std::endl;
    }

    return u;
}