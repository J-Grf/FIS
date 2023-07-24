#include <iostream>
#include <string>
#include <vector>

#include "matrix.hpp"
#include "EVMethods.hpp"

int main (int argc, char* argv[]) {

    if(argc != 2) {
        std::cerr << "Provide method and path to input-File" << std::endl;
        return -1;
    }

    Matrix A{argv[2]};
    readFromFile(A);

    // initial guess
    std::vector<double> x(1/sqrt(A.getDim()), A.getDim());
    if(!static_cast<bool>(strcmp(argv[1], "POW"))) {
        

    } else if (!static_cast<bool>(strcmp(argv[1], "LANC"))) {


    } else {
        std::cerr << "Provide method 'POW' (Power iteration) or 'LANC' (Lanczos method)! " << std::endl;
        return -1;
    }
    return 0;
}