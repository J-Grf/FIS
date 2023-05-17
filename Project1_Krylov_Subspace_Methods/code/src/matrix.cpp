#include <cassert>
#include "matrix.hpp"

Matrix::Matrix (const string _inputDir ) : inputDir(_inputDir) {}

void Matrix::print() const {

    std::string sym = "";
    if(sym_flag == 's') {
        sym = "symmetric";
    } else if(sym_flag == 'n') {
        sym = "non-symmetric";
    } else {
        sym = "undefined TYPE";
    }
    
    std::string JMS, VMS;
    for(size_t i = 0; i < array_size; i++) {
        JMS += std::to_string(JM[i]) + "\n";
        VMS += std::to_string(VM[i]) + "\n";
    } 

    std::cout << "Matrix with input: " << inputDir << std::endl
              << "has type: " << sym << std::endl
              << "has dimesions: " << dim << std::endl
              << "array_size: " << array_size << std::endl
              << "----------------JM-array:-------" << std::endl
              << JMS << std::endl
              << "----------------VM-array:-------" << std::endl
              << VMS << std::endl;

}

void readFromFile (Matrix& A) {
    std::ifstream in(A.inputDir);
    
    if(in.is_open()) {
        in >> A.sym_flag; assert((A.sym_flag == 's' || A.sym_flag == 'n') && "undefined matrix-TYPE");
        in >> A.dim;
        in >> A.array_size;
        A.JM.resize(A.array_size);
        A.VM.resize(A.array_size);
        for(size_t i = 0; i < A.array_size; i++) {
            in >> A.JM[i] >> A.VM[i];
        }
    }
    in.close();
}

std::vector<double> vectorProductCSR (const Matrix& A, const std::vector<double> b) {

}

std::vector<double> vectorProductCSC (const Matrix& A, const std::vector<double> b) {

}

std::vector<double> vectorProduct (const Matrix& A, const std::vector<double> b) {

}