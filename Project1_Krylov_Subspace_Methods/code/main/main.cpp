#include <iostream>
#include <cassert>

#include "matrix.hpp"

int main (int argc, char *argv[]) {

    Matrix test("../../material/test.txt");
    readFromFile(test);
    test.print();

    Matrix gmres("../../material/gmres_test_msr.txt");
    Matrix cg("../../material/cg_test_msr.txt");

    readFromFile(gmres);
    readFromFile(cg);
    
    // debug
    // gmres.print();
}