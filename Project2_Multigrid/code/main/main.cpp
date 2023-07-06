#include <iostream>
#include <string>
#include <cmath>
#include <chrono>
#include <iomanip>

#include "multigrid.hpp"

int main (int argc, char *argv[]) {
    
    if(argc == 7) {
        std::cout << argv[1] << std::endl;
        const std::string strArg1 = std::string(argv[1]);
        const char Arg2 = argv[2][0];
        if(strArg1 != "sg" && strArg1 != "mg" && strArg1 != "test"){
            std::cerr << "method unkown, specify either single- or multigrid or test to test Prolongation and Restriciton! " << std::endl;
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

        bool gridIO = false;
#ifndef DISABLEIO
        gridIO = true;
#endif
        
        const std::vector<std::pair<double, double>> grid = getGrid(N, gridIO, "grid.txt");

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

#ifndef DISABLEIO
            printSolution(u_ex, "u_ex.txt");
#endif
            MG mg(N, n, false);
            mg.setStaticVariables(gamma,nu1,nu2);

            using namespace std::chrono;
            
            high_resolution_clock::time_point t1 = high_resolution_clock::now();
            mg.MG_Algorithm(n, u0, mg.f_rhs);
            high_resolution_clock::time_point t2 = high_resolution_clock::now();
            auto time = duration_cast<nanoseconds>(t2 - t1);

            printf("MG_Algorithm required:  %.5f seconds.\n", time.count() * 1e-9);

            std::ofstream out;
            out.open("../../rawData/timings.txt", std::ios::app);
            /* out << std::left << std::setw(3) << "n" << std::setw(17) << "time in seconds" 
                << std::setw(5) << "nu1" << std::setw(5) << "nu2" << std::setw(7) << "gamma" << std::endl; */
            out << std::left << std::setw(3) << n << std::setw(17) << time.count() * 1e-9 
                << std::setw(5) << nu1 << std::setw(5) << nu2 << std::setw(7) << gamma << std::endl;
            out.close();

#ifndef DISABLEIO
            printSolution(u0, "u.txt");
#endif

            double maxError = -1;
            for( size_t i = 1; i < N; i++) {
                for(size_t j = 1; j < N; j++) {
                    double tmp = abs(u0[i][j] - u_ex[i][j]);
                    if(tmp > maxError)
                        maxError = tmp;
                }
            }
            std::cout << "maxError MG: " << maxError << std::endl;
        
        } else if (!static_cast<bool>(strcmp(argv[1], "test"))) {

            std::cout << "Testing Restriction" << std::endl;
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
            std::cout << "maxError Prolongation: " << maxError << std::endl;
        } 
        
    } else {
        std::cerr << "provide method <sg, mg, test>, choose <'N' N,'n' n> grid parameter (N = 2^n) and iteration parameters <nu1> <nu2> <gamma> " << std::endl;
        return -1; 
    }
    
    
    return 0;
}