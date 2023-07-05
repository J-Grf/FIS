#include <iostream>
#include <string>
#include <cmath>

#include "multigrid.hpp"

int main (int argc, char *argv[]) {
    
    switch (argc) {
        case 1:
            std::cerr << "provide method, choose 'N' or 'n' and specifiy grid parameter N = 2^n and max It nu1 and nu2 + gamma -> cycle type!" << std::endl;
            return -1;
        case 7: {
            std::cout << argv[1] << std::endl;
            const std::string strArg1 = std::string(argv[1]);
            const char Arg2 = argv[2][0];
            if(strArg1 != "sg" && strArg1 != "mg"){
                std::cerr << "method unkown, specify either single- or multigrid! " << std::endl;
                return -1;
            }

            size_t n = 0;
            double N = 0.0; 
            if(Arg2 == 'N') {
                N = std::stod(argv[3]);
            } else if(Arg2 == 'n') {
                n = static_cast<size_t>(std::stoi(argv[3]));
                N = pow(2, n);
            } else {
                std::cerr << "Option " << Arg2 << " unkown, specify either n or N for N = 2^n!" << std::endl;
                return -1; 
            }
            std::cout << "grid points: " << N+1 << std::endl;
            
            const size_t vecSize = static_cast<size_t>(N) + 1;
            const size_t nu1 = static_cast<size_t>(std::stoi(argv[4]));
            const size_t nu2 = static_cast<size_t>(std::stoi(argv[5]));
            const size_t gamma = static_cast<size_t>(std::stoi(argv[6]));  
            m_type u0(vecSize, std::vector<double>(vecSize, 0.0));
            m_type u_ex(vecSize, std::vector<double>(vecSize, 0.0));

            const std::vector<std::pair<double, double>> grid = getGrid(N, true, "grid.txt");

            if(!static_cast<bool>(strcmp(argv[1], "sg"))) {
                
                //Loop over inner points
                m_type f(vecSize, std::vector<double>(vecSize, 0.0));
                double tmp; 
                double initError = 0.0;
                for( size_t i = 1; i < N; i++) {
                    for(size_t j = 1; j < N; j++) {
                        tmp = sin(2 * M_PI * grid[i * vecSize + j].first) * sin(2 * M_PI * grid[i * vecSize + j].second);
                        f[i][j] = 8 * pow(M_PI, 2) * tmp;
                        u_ex[i][j] = tmp;
                        if(abs(u_ex[i][j]) > initError)
                            initError = abs(u_ex[i][j]);
                    }
                }
                printSolution(u_ex, "u_ex.txt");

                m_type u = GaussSeidel(u0, u_ex, f, nu1, N);

                printSolution(u, "u_est.txt");

                std::cout << "initial error: " << initError << std::endl;

            } else if (!static_cast<bool>(strcmp(argv[1], "mg"))) {
                
                u_ex = getExactSolution(N, grid);
                printSolution(u_ex, "u_ex.txt");
                MG mg(N, n, false);
                mg.setStaticVariables(gamma,nu1,nu2);
                mg.MG_Algorithm(n, u0, mg.f_rhs);
                printSolution(u0, "u.txt");

                double maxError = -1;
                for( size_t i = 1; i < N; i++) {
                    for(size_t j = 1; j < N; j++) {
                        double tmp = abs(u0[i][j] - u_ex[i][j]);
                        if(tmp > maxError)
                            maxError = tmp;
                    }
                }
                std::cout << "maxError MG: " << maxError << std::endl;

                /* std::cout << "Testing Restriction" << std::endl;
                //Test restriction
                MG mg(N, n, true);
                double N_c = N/2;
                std::cout << "N_c grid points: " << N_c << std::endl;
                const std::vector<std::pair<double, double>> grid_c = getGrid(N_c, true, "grid_coarse.txt");
                
                m_type u_ex = getExactSolution(N, grid);
                printSolution(u_ex, "u_ex.txt");
                m_type u_coarse = mg.Restriction(u_ex, N_c);
                
                double maxError = -1;
                for( size_t i = 1; i < N_c; i++) {
                    for(size_t j = 1; j < N_c; j++) {
                        double tmp = abs(u_coarse[i][j] - u_ex[2*i][2*j]);
                        if(tmp > maxError) {
                            //std::cout << "at: " << i << " " << j << "tmp: " << tmp << std::endl;
                            maxError = tmp;
                        }
                    }
                }
                std::cout << "maxError Restriction: " << maxError << std::endl;

                std::cout << "Test Prolongation" << std::endl;
                double N_f = 2 * N;
                const std::vector<std::pair<double, double>> grid_f = getGrid(N_f, true, "grid_fine.txt");
                m_type u_ex_f = getExactSolution(N_f, grid_f);
                printSolution(u_ex_f, "u_ex_fine.txt");
                m_type u_fine = mg.Prolongation(u_ex, N);
                maxError = -1;
                for( size_t i = 1; i < N_f; i++) {
                    for(size_t j = 1; j < N_f; j++) {
                        double tmp = abs(u_fine[i][j] - u_ex_f[i][j]);
                        if(tmp > maxError)
                            maxError = tmp;
                    }
                }
                std::cout << "maxError Prolongation: " << maxError << std::endl; */
            } 
            break;
        }
        default:
            std::cerr << "unkown number of arguments, just provide method (sg/mg), choose 'N' or 'n' and specifiy grid parameter N = 2^n  and max It nu!" << std::endl;
            return -1;
    }
    
    return 0;
}