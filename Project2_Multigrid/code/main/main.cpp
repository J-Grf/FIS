#include <iostream>
#include <string>
#include <cmath>

#include "grid.hpp"
#include "gaussseidel.hpp"

int main (int argc, char *argv[]) {
    
    switch (argc) {
        case 1:
            std::cerr << "provide method, grid parameter n -> N = 2^n and the max. number of iterations !" << std::endl;
            return -1;
        case 4: {
            std::cout << argv[1] << std::endl;
            const std::string strArg = std::string(argv[1]);
            if(strArg != "sg" && strArg != "mg"){
                std::cerr << "method unkown, specify either single- or multigrid! " << std::endl;
                return -1;
            }
            // const size_t n = static_cast<size_t>(std::stoi(argv[2]));
            const double N = std::stod(argv[2]); // pow(n,2);
            const size_t nu = static_cast<size_t>(std::stoi(argv[3]));

            m_type u0(N, std::vector<double>(N, 0.0));
            m_type f(N, std::vector<double>(N, 0.0)), u_ex(N, std::vector<double>(N, 0.0));
            const std::vector<std::pair<double, double>> grid = getGrid(N);
            
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

            if(argv[1] == "sg") {

            } else if (argv[2] == "mg") {

            } 
            break;
        }
        default:
            std::cerr << "unkown number of arguments, just provide method (sg/mg), grid parameter n -> N = 2^n and the max. number of iterations !" << std::endl;
            return -1;
    }
    
    return 0;
}