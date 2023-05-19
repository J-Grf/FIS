#ifndef _MATRIX_INCLUDE_
#define _MATRIX_INCLUDE_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

class Matrix {
    char sym_flag = -1;
    size_t dim = 0;
    size_t array_size = 0;
    
    std::vector<int> JM = {};
    std::vector<double> VM = {};

    const std::string inputDir;

public:
    Matrix(const std::string _inputDir);
    size_t getDim() const { return dim; };

    void print() const;

friend void readFromFile (Matrix&);
friend std::vector<double> vectorProduct (const Matrix&, const std::vector<double>&);
};

void readFromFile(Matrix& A);

//vectorProducts in different formats
std::vector<double> vectorProductCSR (const size_t dim, const std::vector<size_t>& IA, const std::vector<size_t>& J,
const std::vector<double>& V, const std::vector<double>& x);

std::vector<double> vectorProductCSC (const size_t dim, const std::vector<size_t>& IA, const std::vector<size_t>& J,
const std::vector<double>& V, const std::vector<double>& x);

std::vector<double> vectorProduct(const Matrix& A, const std::vector<double>& x);

//helper functions
inline double dotP(const std::vector<double> a, const std::vector<double> b);
inline double norm2(const std::vector<double> a);

#endif