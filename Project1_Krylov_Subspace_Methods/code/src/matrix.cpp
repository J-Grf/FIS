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

std::vector<double> vectorProductCSR (const size_t dim, const std::vector<double>& IA, const std::vector<double>& J,
const std::vector<double>& V, const std::vector<double> x) {
    assert(dim == x.size() && "dimensions of matrix A and vector x differ!");
    
    std::vector<double> y(dim, 0.0);
    
    //TO CSR-format
    

    //for(size_t i = 0; i < A.dim + 1; i++) {}
    // IA = A.JM[0:A.dim]
    // J = A.JM[A.dim+1 : A.array_size]
    // V = A.VM[A.dim+1 : A.array_size]
    
    int i1 = -1;
    int i2 = -1;
    std::vector<double> tmpV;
    std::vector<double> tmpX;
    int tmplength = 0;
    size_t idx = 0;
    for(size_t i = 0; i < dim; i++) {
        
        i1 = IA[i];
        i2 = IA[i+1]-1;

        tmplength = i2 - i1;
        tmplength++;
        
        //std::cout << "i1: " << i1 << " i2: " << i2 << " tmplength : " << tmplength << std::endl;
        
        if (tmplength == 0 ) 
            continue;
        
        tmpV.resize(tmplength);
        tmpX.resize(tmplength);
        
        for(size_t j = 0; j < tmplength; j++) {
            //idx = A.dim + 1 + i1;
            //std::cout << "idx: " << idx << std::endl;
            tmpV[j] = V[i1];
            //std::cout << "i1: " << i1 << std::endl;
            //std::cout << "tmpV[" << j << "]: " << tmpV[j] << std::endl;
            tmpX[j] = x[J[i1]]; 
            //std::cout << "J[" << i1 << "]: " << J[i1] << std::endl;
            //std::cout << " tmpX[ " << j << "]: " << tmpX[j] << std::endl;
            i1++;
            if( i1 == J.size())
                break;
        }
        
        y[i] = dot(tmpV, tmpX);
        tmpV.clear();
        tmpX.clear();
    }

    return y;
}

std::vector<double> vectorProductCSC (const size_t dim, const std::vector<double>& IA, const std::vector<double>& J,
const std::vector<double>& V, const std::vector<double> x) {
    assert(dim == x.size() && "dimensions of matrix A and vector x differ!");
    
    std::vector<double> y(dim, 0.0);

    int i1 = -1;
    int i2 = -1;
    int tmplength = 0;
    size_t idx = 0;
    for(size_t i = 0; i < dim; i++) {
        
        i1 = IA[i];
        i2 = IA[i+1]-1;

        tmplength = i2 - i1;
        tmplength++;
        
        //std::cout << "i1: " << i1 << " i2: " << i2 << " tmplength : " << tmplength << std::endl;
        
        if (tmplength == 0 ) 
            continue;
        
        for(size_t j = 0; j < tmplength; j++) {
            y[J[i1]] += V[i1] * x[i];
            i1++;
            if( i1 == J.size())
                break;
        }
    }

    return y;
}

std::vector<double> vectorProduct (const Matrix& A, const std::vector<double> x) {
    assert(A.dim == x.size() && "dimensions of matrix A and vector x differ!");

    std::vector<double> y(A.dim, 0);

    // 1) compute contribution of diagonal elements
    for(size_t i = 0; i < A.dim; i++) {
        y[i] = A.VM[i] * x[i];
    }


    // 2) compute contribution of off-diagonal elements
    if(A.sym_flag == 's') {
        // Algorithm for symmetric matrix in the form of D * x + U^T * x + U * x = y
        
        // CSR for U

        // CSC for U^T
    } else {
        //CSR for off-diagonal entries
        //Create CSR-like format for off-diagonal entries

        std::vector<double> IA;
        std::vector<double> J;
        std::vector<double> V;
        int shift = A.dim + 2;
        for(size_t i = 0; i < A.dim + 1; i++) {

            IA.push_back(A.JM[i] - shift);
        }

        for(size_t i = A.dim + 1; i < A.array_size; i++) {
            J.push_back(A.JM[i]);
            V.push_back(A.VM[i]);
        }
        std::vector<double> yoff = vectorProductCSR(A.dim, IA, J, V, x);
        
        for(size_t i = 0; i < A.dim; i++) {
            y[i] += yoff[i];
        }
    }

    return y;
}

double dot(const std::vector<double> a, const std::vector<double> b){
    assert(b.size()==a.size());
    double res = 0;
    for (size_t i = 0; i < a.size(); i++) {
        res += a[i]*b[i];
    }
    return res;
}