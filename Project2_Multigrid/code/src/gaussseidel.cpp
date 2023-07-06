#include "gaussseidel.hpp"
#include <iostream>
#include <fstream>

m_type GaussSeidel(m_type& u, const m_type& u_ex, const m_type& f, const size_t nu, const size_t N) {
    std::ofstream out, out2;
    out.open("maxError.txt");
    out2.open("infNorm.txt");

    const double hsq = pow(1 / static_cast<double>(N) ,2);
    /*
        - u[i-1][j] and u[i][j-1] at k
        - u[i+1][j] and u[i][j+1] at k - 1 from previous iteration!
    */
    m_type u_old = u;
    size_t k = 0;
    double infNorm;
    double max_error;
    do {
        k++;
        double diffU, diffUex;
        max_error = -1;
        infNorm = -1;
        for(size_t i = 1; i < N; i++) {
            for(size_t j = 1; j < N; j++) {
                u[i][j] = 0.25 * (hsq * f[i][j] + u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1]);
                
                diffU = abs(u[i][j] - u_old[i][j]);
                diffUex = abs(u_ex[i][j] - u[i][j]);
                if(diffU > infNorm) 
                    infNorm = diffU;
                if(diffUex > max_error)
                    max_error = diffUex;
            }
        }
        u_old = u;
        out << max_error << std::endl;
        out2 << infNorm << std::endl;
        std::cout << k << "  " << infNorm << "  " << max_error << std::endl;
        
        if (k == nu)
            break;
    } while(infNorm > eps);
    
    std::cout << "maximum converged error: " << max_error << std::endl;
    out.close();
    out2.close();

    return u;
}

m_type GaussSeidel(m_type& u, const m_type& f, const size_t nu, const size_t N) {
    assert(f.size() == u.size());
    const double hsq = pow(1 / static_cast<double>(N) ,2);

    /*
        - u[i-1][j] and u[i][j-1] at k
        - u[i+1][j] and u[i][j+1] at k - 1 from previous iteration!
    */
    
    for(size_t k = 0; k < nu; k++) {
        //excluding boundary nodes
        for(size_t i = 1; i < N; i++) {
            for(size_t j = 1; j < N; j++) {
                u[i][j] = 0.25 * (hsq * f[i][j] + u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1]);
            }
        }
    }
    
    return u;
}