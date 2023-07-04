#include <iostream>
#include <string>
#include <cmath>

#include "multigrid.hpp"

int main (int argc, char *argv[]) {
    
    switch (argc) {
        case 1:
            std::cerr << "provide method, choose 'N' or 'n' and specifiy grid parameter N = 2^n and max It nu!" << std::endl;
            return -1;
        case 5: {
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
                std::cout << "N grid points: " << N << std::endl;
            } else {
                std::cerr << "Option " << Arg2 << " unkown, specify either n or N for N = 2^n!" << std::endl;
                return -1; 
            }

            size_t nu = static_cast<size_t>(std::stoi(argv[4]));  
            m_type u0(N, std::vector<double>(N, 0.0));
            m_type u_ex(N, std::vector<double>(N, 0.0));
            const std::vector<std::pair<double, double>> grid = getGrid(N, true, "grid.txt");
            //Loop over inner points
            m_type f(N, std::vector<double>(N, 0.0));

            if(!static_cast<bool>(strcmp(argv[1], "sg"))) {
            
                double tmp;
                double initError = -1;
                for( size_t i = 1; i < N - 1; i++) {
                    for(size_t j = 1; j < N - 1; j++) {
                        tmp = sin(2 * M_PI * grid[i * N + j].first) * sin(2 * M_PI * grid[i * N + j].second);
                        f[i][j] = 8 * pow(M_PI, 2) * tmp;
                        u_ex[i][j] = tmp;
                        if(abs(u_ex[i][j]) > initError)
                            initError = abs(u_ex[i][j]);
                    }
                }

                m_type u = GaussSeidel(u0, u_ex, f, nu);

                std::cout << "initial error: " << initError << std::endl;

            } else if (!static_cast<bool>(strcmp(argv[1], "mg"))) {
                std::cout << "Testing Restriction" << std::endl;
                //Test restriction
                MG mg(N, n);
                double N_c = N/2;
                std::cout << "N_c grid points: " << N_c << std::endl;
                const std::vector<std::pair<double, double>> grid_c = getGrid(N_c, true, "grid_coarse.txt");
                
                m_type u_ex = getExactSolution(N, grid);
                printExactSolution(u_ex, "u_ex.txt");
                m_type u_coarse = mg.Restriction(u_ex, N_c);
                
                double maxError = -1;
                for( size_t i = 1; i < N_c - 1; i++) {
                    for(size_t j = 1; j < N_c - 1; j++) {
                        double tmp = abs(u_coarse[i][j] - u_ex[2*i][2*j]);
                        if(tmp > maxError) {
                            std::cout << "at: " << i << " " << j << "tmp: " << tmp << std::endl;
                            maxError = tmp;
                        }
                    }
                }
                std::cout << "maxError Restriction: " << maxError << std::endl;

                std::cout << "Test Prolongation" << std::endl;
                double N_f = 2 * N;
                const std::vector<std::pair<double, double>> grid_f = getGrid(N_f, true, "grid_fine.txt");
                m_type u_ex_f = getExactSolution(N_f, grid_f);
                printExactSolution(u_ex_f, "u_ex_fine.txt");
                m_type u_fine = mg.Prolongation(u_ex, N);
                maxError = -1;
                for( size_t i = 1; i < N_f - 1; i++) {
                    for(size_t j = 1; j < N_f - 1; j++) {
                        double tmp = abs(u_fine[i][j] - u_ex_f[i][j]);
                        if(tmp > maxError)
                            maxError = tmp;
                    }
                }
                std::cout << "maxError Prolongation: " << maxError << std::endl;

            } 
            break;
        }
        default:
            std::cerr << "unkown number of arguments, just provide method (sg/mg), choose 'N' or 'n' and specifiy grid parameter N = 2^n  and max It nu!" << std::endl;
            return -1;
    }
    
    return 0;
}