#include <iostream>

#include "matrix.hpp"
#include "gmres.hpp"
#include "cg.hpp"

int main (int argc, char *argv[]) {
    //Read gmres-test matrix
    if(argc == 1) {
        std::cerr << "provide arguments! " << std::endl;
        return -1;
    }
        
    if( std::string(argv[1]) == "GMRES" ) {
        if(argc != 5) {
            std::cerr << "Provide 3 more arguments when using GMRES: <path> <restartParameter> <preconditioner> " << std::endl;
            return -1;
        }

        Matrix gmresTest(argv[2]);
        readFromFile(gmresTest);

        //Prescribe solution vector and determine right-hand-side
        const std::vector<double> x(gmresTest.getDim(), 1.0);
        const std::vector<double> b = vectorProduct(gmresTest, x);
        const std::vector<double> x0(gmresTest.getDim());
        const size_t m = static_cast<size_t>(std::stoi(argv[3]));

        const std::vector<double> res = GMRES_Res(gmresTest, x0, b, m, StringToPre.at(argv[4]));
        
        /* for(size_t i = 0; i < res.size(); i++){
            std::cout << "x[" << i << "]: " << res[i] << std::endl;
        } */
    } else if ( std::string(argv[1]) == "CG" ) {
        if(argc != 4) {
            std::cerr << "Provide 2 more argument when using CG: <path> <iterations>" << std::endl;
            return -1;
        }
        Matrix cgTest(argv[2]);
        readFromFile(cgTest);

        //Prescribe solution vector and determine right-hand-side
        const std::vector<double> x(cgTest.getDim(), 1.0);
        const std::vector<double> b = vectorProduct(cgTest, x);
        const std::vector<double> x0(cgTest.getDim());
        const size_t m = static_cast<size_t>(std::stoi(argv[3]));

        auto res = CGMethod(cgTest, b, x0, m);
        
        /*for(size_t i = 0; i < res.first.size(); i++) {
            std::cout << "x[" << i << "]: " << res.first[i] << std::endl;
        }*/
    } else {
        std::cerr << "Input argument is unkown, choose between 'GMRES' or 'CG' " << std::endl;
        return -1;
    }

    return 0;
}