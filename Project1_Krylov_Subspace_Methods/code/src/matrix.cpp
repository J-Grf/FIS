#include <cassert>
#include <cmath>
#include "preconditioner.hpp"
#include "matrix.hpp"

Matrix::Matrix (const std::string& _inputDir ) : inputDir(_inputDir) {}

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
    JMS += "[";
    VMS = JMS;
    for(size_t i = 0; i < array_size; i++) {
        JMS += std::to_string(JM[i]) + ", ";
        VMS += std::to_string(VM[i]) + ", ";
    } 
    JMS += "]";
    VMS += "]";

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
    std::cout << "Reading from file: " << A.inputDir << std::endl;
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

std::vector<double> vectorProductCSR (const size_t dim, const std::vector<size_t>& IA, const std::vector<size_t>& J,
const std::vector<double>& V, const std::vector<double>& x) {
    assert(dim == x.size() && "dimensions of matrix A and vector x differ!");
    
    std::vector<double> y(dim, 0.0);
    
    int i1 = -1;
    int i2 = -1;
    std::vector<double> tmpV;
    std::vector<double> tmpX;
    int tmplength = 0;
    for(size_t i = 0; i < dim; i++) {
        
        i1 = IA[i];
        i2 = IA[i+1]-1;

        tmplength = i2 - i1;
        tmplength++;
        
        if (tmplength <= 0 ) 
            continue;
        
        tmpV.resize(tmplength);
        tmpX.resize(tmplength);
        
        for(size_t j = 0; j < tmplength; j++) {
            if( i1 == J.size())
                break;
            tmpV[j] = V[i1];
            tmpX[j] = x[J[i1]]; 
            i1++;
        }
        
        y[i] = dotP(tmpV, tmpX);
        tmpV.clear();
        tmpX.clear();
    }

    return y;
}

std::vector<double> vectorProductCSC (const size_t dim, const std::vector<size_t>& IA, const std::vector<size_t>& J,
const std::vector<double>& V, const std::vector<double>& x) {
    assert(dim == x.size() && "dimensions of matrix A and vector x differ!");
    
    std::vector<double> y(dim, 0.0);

    int i1 = -1;
    int i2 = -1;
    int tmplength = 0;
    for(size_t i = 0; i < dim; i++) {
        
        i1 = IA[i];
        i2 = IA[i+1]-1;

        tmplength = i2 - i1;
        tmplength++;
        
        if (tmplength <= 0 ) 
            continue;
        
        for(size_t j = 0; j < tmplength; j++) {
            if( i1 == J.size())
                break;
            y[J[i1]] += V[i1] * x[i];
            i1++;
        }
    }

    return y;
}

std::vector<double> vectorProduct (const Matrix& A, const std::vector<double>& x) {
    assert(A.dim == x.size() && "dimensions of matrix A and vector x differ!");

    std::vector<double> y(A.dim, 0);

    // 1) compute contribution of diagonal elements
    for(size_t i = 0; i < A.dim; i++) {
        y[i] = A.VM[i] * x[i];
    }

    //2) Create CSR-like format for off-diagonal entries
    std::vector<double> V;
    std::vector<size_t> IA, J;
    const int shift = A.dim + 2;
    for(size_t i = 0; i < A.dim + 1; i++) {

        IA.push_back(A.JM[i] - shift);
    }

    for(size_t i = A.dim + 1; i < A.array_size; i++) {
        J.push_back(A.JM[i] - 1);
        V.push_back(A.VM[i]);
    }

    // 2) compute contribution of off-diagonal elements
    if(A.sym_flag == 's') {
        // Algorithm for symmetric matrix in the form of D * x + U^T * x + U * x = y
        // CSR for U
        std::vector<double> yoffUpper = vectorProductCSR(A.dim, IA, J, V, x);
        // CSC for U^T
        std::vector<double> yoffLower = vectorProductCSC(A.dim, IA, J, V, x);

        for(size_t i = 0; i < A.dim; i++) {
            y[i] += yoffUpper[i] + yoffLower[i];
        }

    } else {
        //CSR for off-diagonal entries
        std::vector<double> yoff = vectorProductCSR(A.dim, IA, J, V, x);
        
        for(size_t i = 0; i < A.dim; i++) {
            y[i] += yoff[i];
        }
    }

    return y;
}

std::vector<double> backwardSubMSR(const Matrix& A, const std::vector<double>& b, const size_t m) {
    assert(A.getDim() == b.size() && "Dimensions of A and b do not coincide, backwardSub not possible!");
    assert(A.sym_flag == 'n' && "backwardSub only applicable for upper triangular Matrix (non-sym)");

    std::vector<double> x(m, 0);
    size_t i1 = m - 1;
    size_t i2 = 0;
    size_t k = 1;
    int tmplength = 0;

    //Get first entry of x, bottom of upper triangle
    x[i1] = b[i1] / A.VM[i1]; //here m should correspond to A.dim

    //next compute all the following rows
    for(int i = m - 2; i >= 0; i--, k++) {
        x[i] = b[i];

        //find start and end index of off-diagonal elements in rows
        i1 = A.JM[A.dim - k - 1] - 1;
        i2 = A.JM[A.dim - k] - 1;
        tmplength = i2 - i1;
        
        if (tmplength <= 0 ) 
            continue;

        for(size_t j = 0; j < tmplength; j++) {
            // extract upper diagonal entries by checking if column index is larger than the row index
            if(A.JM[i1] - 1 > i) {
                //get off-diagonal element and multiply with corresponding x element at column index of A
                x[i] -= A.VM[i1] * x[A.JM[i1] - 1];
            }
            i1++;
        }
        x[i] /= A.VM[i];
    }
    return x;
}

std::vector<double> forwardSubMSR(const Matrix& A, const std::vector<double>& b, const size_t m) {
    assert(A.getDim() == b.size() && "Dimensions of A and b do not coincide, forwardSub not possible!");
    assert(A.sym_flag == 'n' && "forwardSub only applicable for lower triangular Matrix (non-sym)");

    std::vector<double> x(m, 0);
    size_t i1 = 0;
    size_t i2 = 0;
    int tmplength = 0;

    //Get first entry of x, top of lower triangle
    x[i1] = b[i1] / A.VM[i1]; //here m should correspond to A.dim

    //next compute all the following rows
    for(size_t i = 1; i < A.dim; i++) {
        x[i] = b[i];

        //find start and end index of off-diagonal elements in rows
        i1 = A.JM[i] - 1;
        i2 = A.JM[i + 1] - 1;
        tmplength = i2 - i1;
        
        if (tmplength <= 0 ) 
            continue;

        for(size_t j = 0; j < tmplength; j++) {
            // extract lower diagonal entries by checking if column index is smaller than the row index
            if(A.JM[i1] - 1 < i) {
                //get off-diagonal element and multiply with corresponding x element at column index of A
                x[i] -= A.VM[i1] * x[A.JM[i1] - 1];
            }
            i1++;
        }
        x[i] /= A.VM[i];
    }
    return x;
}