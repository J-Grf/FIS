#include <iostream>
#include <cassert>
#include <chrono>

#include "matrix.hpp"
#include "gmres.hpp"
#include "cg.hpp"

int main (int argc, char *argv[]) {

    /* Matrix test("../../material/test3.txt");
    readFromFile(test);

    ILUout A = Ilu(test);
    for(size_t i = 0; i < A.PV.size(); i++) {
        std::cout << "PV[" << i << "]: " << A.PV[i] << std::endl;
    }

    for(size_t i = 0; i < A.PJM.size(); i++) {
        std::cout << "PJM[" << i << "]: " << A.PJM[i] << std::endl;
    }

    for(size_t i = 0; i < A.JD.size(); i++) {
        std::cout << "PJD[" << i << "]: " << A.JD[i] << std::endl;
    }

    exit(0); */

    //Read gmres-test matrix
    if(argc == 1) {
        std::cerr << "provide arguments! " << std::endl;
        return -1;
    }
        
    if( std::string(argv[1]) == "GMRES" ) {
        if(argc != 5) {
            std::cerr << "Provide 2 more arguments when using GMRES: <path> <restartParameter> <preconditioner> " << std::endl;
            return -1;
        }

        Matrix gmresTest(argv[2]);
        readFromFile(gmresTest);
        //gmresTest.print();

        //Prescribe solution vector and determine right-hand-side
        const std::vector<double> x(gmresTest.getDim(), 1.0);
        const std::vector<double> b = vectorProduct(gmresTest, x);
        const std::vector<double> x0(gmresTest.getDim());
        const size_t m = static_cast<size_t>(std::stoi(argv[3]));

        auto begin = std::chrono::high_resolution_clock::now();

        const std::vector<double> res = GMRES_Res(gmresTest, x0, b, m, StringToPre.at(argv[4]));
        auto end = std::chrono::high_resolution_clock::now();
        
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf("Time measured: %.3f seconds.\n", elapsed.count() * 1e-9);
        
        for(size_t i = 0; i < res.size(); i++){
            std::cout << "x[" << i << "]: " << res[i] << std::endl;
        }
    } else if ( std::string(argv[1]) == "CG" ) {
        if(argc != 4) {
            std::cerr << "Provide 2 more argument when using CG: <path> <iterations>" << std::endl;
            return -1;
        }
        Matrix cgTest(argv[2]);
        readFromFile(cgTest);
        //cgTest.print();

        //Prescribe solution vector and determine right-hand-side
        const std::vector<double> x(cgTest.getDim(), 1.0);
        const std::vector<double> b = vectorProduct(cgTest, x);
        const std::vector<double> x0(cgTest.getDim());
        const size_t m = static_cast<size_t>(std::stoi(argv[3]));

        auto begin = std::chrono::high_resolution_clock::now();
        auto res = CGMethod(cgTest, b, x0, m);
        auto end = std::chrono::high_resolution_clock::now();
        
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf("Time measured: %.3f seconds.\n", elapsed.count() * 1e-9);
        
        for(size_t i = 0; i < res.first.size(); i++) {
            std::cout << "x[" << i << "]: " << res.first[i] << std::endl;
        }
    } else {
        std::cerr << "Input argument is unkown, choose between 'GMRES' or 'CG' " << std::endl;
        return -1;
    }

    return 0;

    //For debugging with test matrix
    /* Matrix test("../../material/test.txt");
    readFromFile(test);
    test.print();
    std::vector<double> t {1.123123, 2.123213, 0.1323, 4.2770123E-13};
    std::vector<size_t> IA1 {0, 2, 4, 5, 7};
    std::vector<size_t> J1 {0, 2, 1, 2, 2, 2, 3};
    std::vector<double> V1 {1, 4, 2, 2, 1, 4, 2};

    Matrix test("../../material/test_sym.txt");
    readFromFile(test);
    test.print();
    std::vector<double> t {1.123123, 2.123213, 0.1323, 4.2770123E-13};
    std::vector<size_t> IA1 {0, 2, 4, 7, 8};
    std::vector<size_t> J1 {0, 2, 1, 2, 0, 1, 2, 3};
    std::vector<double> V1 {1, 4, 2, 2, 4, 2, 1, 2};

    std::vector<double> b = vectorProductCSR(test.getDim(), IA1, J1, V1, t);
    
    std::vector<size_t> IA2 {0, 1, 2, 6, 7};
    std::vector<size_t> J2 {0, 1, 0, 1, 2, 3, 3};
    std::vector<double> V2 {1, 2, 4, 2, 1, 4, 2}; 
    
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

    */

    /* Matrix gmres("../../material/gmres_test_msr.txt");
    Matrix cg("../../material/cg_test_msr.txt");

    readFromFile(gmres);
    readFromFile(cg); */
    
    // debug
    // gmres.print();

    //debug
    //construct orthonormal basis
    /* Matrix test("../../material/test3.txt");
    readFromFile(test);
    test.print();

    std::vector<double> xp {1, 1, 1};
    std::vector<double> b = vectorProduct(test, xp);
    std::vector<double> r0 = b;
    auto res = gramSchmidt(test, r0, 2);
    printGM(res, 2);
    std::cout << "orthogonality check: " << std::endl << std::boolalpha
              << checkOrthogonality(res.first.at(0), res.first.at(1)) << " "
              << checkOrthogonality(res.first[1], res.first[2]) << " "
              << checkOrthogonality(res.first[0], res.first[2]) << " " << std::endl; 

    //initial guess
    std::vector<double> x0 {0, 0 ,0};

    //std::vector<double> t {1.123123, 2.123213, 0.1323};
    std::vector<double> x = MR_method(test, b, x0);
    for(size_t i = 0; i < x.size(); i++){
        std::cout << "x[" << i << "]: " << x[i] << std::endl;
    }

    std::vector<double> xG = GMRES_Res(test, x0, b, 3, NONE);
    std::cout << "--------- FINAL -------" << std::endl;
    for(size_t i = 0; i < xG.size(); i++){
        std::cout << "xG[" << i << "]: " << xG[i] << std::endl;
    } */


    //Matrix testlow("../../material/test_lower.txt");
    //readFromFile(testlow);
    //testlow.print();
    //Test backward substitution
    /* matrixType<double> A {{1, 2.3, 0.782, 212}, {0, 1, 0, 7}, {0, 0, 3 , 1}, {0 , 0 , 0 ,4}};
    for(size_t i = 0; i < A.size(); i++) {
        for(size_t j = 0; j < A[0].size(); j++){
            std::cout << "A[" << i << "][" << j << "]: " << A[i][j] << std::endl;
        }
    }
    std::vector<double> b {1.59934, -0.125656, 5.234, 2.13123, 0.0};

    std::vector<double> x = backwardSub(A, b, 4);
    for(size_t i = 0; i < x.size(); i++){
        std::cout << "x[" << i << "]: " << x[i] << std::endl;
    }

    std::vector<double> b_test = VP(A,x);
    for(size_t i = 0; i < b_test.size(); i++){
        std::cout << "b_test[" << i << "]: " << b_test[i] << std::endl;
    } */

    //std::vector<double> r0 {1, 1, 1};
    
    //compute rhs b from prescribed solution vector
    //std::vector<double> xp {1.567, 4.562, 8.312};
    //std::vector<double> b = {43.9, 86.6, 121.4, 147.2};//vectorProduct(testlow, xp);

    
    /* for(size_t i = 0; i < b.size(); i++){
        std::cout << "b[" << i << "]: " << b[i] << std::endl;
    } 

    //auto resCG = CGMethod(test, b, x0, m);
    //check forward and backward substitutions:
    std::vector<double> b1 {1.3,  12.5, 33.1, 193.1};
    std::vector<double> w = forwardSubMSR(test, b1, test.getDim());
    for(size_t i = 0; i < w.size(); i++){
        std::cout << "w[" << i << "]: " << w[i] << std::endl;
    }

    std::vector<double> b2 = {43.9, 86.6, 121.4, 147.2};
    std::vector<double> x = backwardSubMSR(test, b2, test.getDim());
    for(size_t i = 0; i < x.size(); i++){
        std::cout << "x[" << i << "]: " << x[i] << std::endl;
    } */
}