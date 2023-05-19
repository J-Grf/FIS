#include <iostream>
#include <cassert>

#include "matrix.hpp"

int main (int argc, char *argv[]) {

    /* Matrix test("../../material/test.txt");
    readFromFile(test);
    test.print();
    std::vector<double> t {1.123123, 2.123213, 0.1323, 4.2770123E-13};
    std::vector<size_t> IA1 {0, 2, 4, 5, 7};
    std::vector<size_t> J1 {0, 2, 1, 2, 2, 2, 3};
    std::vector<double> V1 {1, 4, 2, 2, 1, 4, 2}; */

    Matrix test("../../material/test_sym.txt");
    readFromFile(test);
    test.print();
    std::vector<double> t {1, 2, 3, 4};
    std::vector<size_t> IA1 {0, 2, 4, 7, 8};
    std::vector<size_t> J1 {0, 2, 1, 2, 0, 1, 2, 3};
    std::vector<double> V1 {1, 4, 2, 2, 4, 2, 1, 2};

    std::vector<double> b = vectorProductCSR(test.getDim(), IA1, J1, V1, t);
    
    /* std::vector<size_t> IA2 {0, 1, 2, 6, 7};
    std::vector<size_t> J2 {0, 1, 0, 1, 2, 3, 3};
    std::vector<double> V2 {1, 2, 4, 2, 1, 4, 2}; */ 
    
    std::vector<double> c = vectorProductCSC(test.getDim(), IA1, J1, V1, t);

    std::vector<double> y = vectorProduct(test, t);
    std::string bS = "[";
    std::string cS = bS;
    std::string yS = bS;
    for(size_t i = 0; i < b.size(); i++) {
        bS += std::to_string(b[i]) + ", "; 
        cS += std::to_string(c[i]) + ", ";
        yS += std::to_string(y[i]) + ", ";
    }
    bS += "]";
    cS += "]";
    yS += "]";
    std::cout << bS << std::endl
              << cS << std::endl
              << yS << std::endl;

    /* Matrix gmres("../../material/gmres_test_msr.txt");
    Matrix cg("../../material/cg_test_msr.txt");

    readFromFile(gmres);
    readFromFile(cg); */
    
    // debug
    // gmres.print();
}