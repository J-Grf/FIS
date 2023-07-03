#include <iostream>
#include <string>
#include <cmath>

#include "grid.hpp"
#include "gaussseidel.hpp"
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
            m_type f(N, std::vector<double>(N, 0.0)), u_ex(N, std::vector<double>(N, 0.0));
            const std::vector<std::pair<double, double>> grid = getGrid(N, true);

            if(!strcmp(argv[1], "sg")) {
                //Loop over inner points
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


                std::cout << "initial error: " << initError << std::endl;
                m_type u = GaussSeidel(u0, u_ex, f, nu);

            } else if (!strcmp(argv[2], "mg")) {
                


            } 
            break;
        }
        default:
            std::cerr << "unkown number of arguments, just provide method (sg/mg), choose 'N' or 'n' and specifiy grid parameter N = 2^n  and max It nu!" << std::endl;
            return -1;
    }
    
    return 0;
}