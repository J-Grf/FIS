#include <iostream>
#include <string>
#include <vector>
#include <chrono>

#include "matrix.hpp"
#include "EVMethods.hpp"

int main (int argc, char* argv[]) {

    if(argc > 4 || argc < 3) {
        std::cerr << "Provide method and path to input-File" << std::endl;
        return -1;
    }

    Matrix A{argv[2]};
    readFromFile(A);
    //A.print();
    using namespace std::chrono;
    // initial guess
    std::vector<double> x(A.getDim(), 1/sqrt(A.getDim()));
    std::ofstream out2;
    if(!static_cast<bool>(strcmp(argv[1], "POW"))) {
        high_resolution_clock::time_point t1 = high_resolution_clock::now();
        const double lambdaMax = powerIteration(A, x, 1E-8);
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        auto time = duration_cast<nanoseconds>(t2 - t1);
        
        std::cout << std::setprecision(17) << std::scientific <<"LambdaMax: " << lambdaMax << std::endl;
        printf("Power Iteration required:  %.5f seconds.\n", time.count() * 1e-9);

    } else if (!static_cast<bool>(strcmp(argv[1], "LANC"))) {
        if(argc != 4) {
            std::cerr << "For Lanczos method, path to input-File and the desired dimensions \
            of the Krylov space have to be specified!" << std::endl;
            return -1;
        }
        const size_t m = static_cast<size_t>(std::stoi(argv[3]));
        high_resolution_clock::time_point t1 = high_resolution_clock::now();
        const double lambdaMax = LanczosMethod(A, x, m);
        high_resolution_clock::time_point t2 = high_resolution_clock::now();
        auto time = duration_cast<nanoseconds>(t2 - t1);

        out2.open("../../rawData/timings_Lanczos.txt", std::ios::app);
        out2 << std::left << std::setw(12) << m << std::setw(12) << time.count() * 1e-9 << std::setw(12) 
             << std::setprecision(17) << std::scientific << lambdaMax << std::endl;
             //abs(lambdaCG - lambdaMax) << std::endl;
        out2.close();
        
        std::cout << std::setprecision(17) << std::scientific << "LambdaMax: " << lambdaMax << std::endl;
        printf("Lanczos Method required:  %.5f seconds.\n", time.count() * 1e-9);

    } else {
        std::cerr << "Provide method 'POW' (Power iteration) or 'LANC' (Lanczos method)! " << std::endl;
        return -1;
    }
    return 0;
}