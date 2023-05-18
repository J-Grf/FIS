#include <iostream>
#include <cassert>

#include "matrix.hpp"

int main (int argc, char *argv[]) {

    Matrix test("../../material/test.txt");
    readFromFile(test);
    test.print();
    std::vector<double> t(test.getDim(), 1.0);
    std::vector<double> IA1 {0, 2, 4, 5, 7};
    std::vector<double> J1 {0, 2, 1, 2, 2, 2, 3};
    std::vector<double> V1 {1, 4, 2, 2, 1, 4, 2};

    std::vector<double> b = vectorProductCSR(test.getDim(), IA1, J1, V1, t);
    
    std::vector<double> IA2 {0, 1, 2, 6, 7};
    std::vector<double> J2 {0, 1, 0, 1, 2, 3, 3};
    std::vector<double> V2 {1, 2, 4, 2, 1, 4, 2}; 
    
    std::vector<double> c = vectorProductCSC(test.getDim(), IA2, J2, V2, t);

    std::vector<double> y = vectorProduct(test, t);
    std::string bS = "[";
    std::string cS = bS;
    std::string yS = bS;
    for(size_t i = 0; i < b.size(); i++) {
        bS += std::to_string(b[i]) + ", "; 
        cS += std::to_string(c[i]) + ", ";
        yS += std::to_string(c[i]) + ", ";
    }
    bS += "]";
    cS += "]";
    yS += "]";
    std::cout << bS << std::endl
              << cS << std::endl
              << yS << std::endl;

    Matrix gmres("../../material/gmres_test_msr.txt");
    Matrix cg("../../material/cg_test_msr.txt");

    readFromFile(gmres);
    readFromFile(cg);
    
    // debug
    // gmres.print();
}