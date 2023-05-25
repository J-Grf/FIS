#ifndef _MATRIX_INCLUDE_
#define _MATRIX_INCLUDE_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
//#include <limits>

template<typename T>
const T eps = std::numeric_limits<T>::epsilon();

const double Eps = 10E-8;
//TODO: make template maybe
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

struct MatrixCoo {
    std::vector<size_t> rows;
    std::vector<size_t> columns;
    std::vector<double> values;

    size_t m = 0;
    size_t n = 0;

    MatrixCoo() = default;
    void detDimensions();
    double& setDiagonal(const size_t index);
    double& at(const size_t i, const size_t j);
    void append(const size_t i, const size_t j, const double value);
    std::vector<double> getDiagonals() const;
    void print() const;
};

//helper functions
inline double dotP(const std::vector<double> a, const std::vector<double> b){
    assert(b.size()==a.size());
    double res = 0;
    for (size_t i = 0; i < a.size(); i++) {
        res += a.at(i)*b.at(i);
    }
    return res;
}

inline double norm2(const std::vector<double> a) {
    return sqrt(dotP(a,a));
}

//conventional vectorProduct (but with transpose)
inline std::vector<double> VP(const std::vector<std::vector<double>>& A, const std::vector<double>& x) {
    assert(A[0].size() == x.size());

    std::vector<double> res(x.size());
    for(size_t i = 0; i < x.size(); i++) {
        for(size_t j = 0; j < x.size(); j++) {
            res[i] += A[j][i] * x[j];
        }
    }
    return res;
}

//implement overloaded operators for std::vector
template<typename T>
std::vector<T>& operator-=(std::vector<T>& v, const std::vector<T>& a) {
    assert(v.size() == a.size());
    for(size_t i = 0; i < v.size(); i++)
        v.at(i) -= a.at(i);
    return v;
}

template<typename T>
std::vector<T>& operator+=(std::vector<T>& v, const std::vector<T>& a) {
    assert(v.size() == a.size());
    for(size_t i = 0; i < v.size(); i++)
        v.at(i) += a.at(i);
    return v;
}

template<typename T>
std::vector<T> operator-(const std::vector<T>& v, const std::vector<T>& a) {
    assert(v.size() == a.size());
    std::vector<T> res(v.size());
    for(size_t i = 0; i < v.size(); i++)
        res[i] = v.at(i) - a.at(i);
    return res;
}

template<typename T>
std::vector<T> operator+(const std::vector<T>& v, const std::vector<T>& a) {
    assert(v.size() == a.size());
    std::vector<T> res(v.size());
    for(size_t i = 0; i < v.size(); i++)
        res[i] = v.at(i) + a.at(i);
    return res;
}

template<typename T>
std::vector<T> operator*(const T& a, const std::vector<T>& v) {
    std::vector<T> res(v.size());
    for(size_t i = 0; i < v.size(); i++) {
        res.at(i) = a * v.at(i);
    }
    return res;
}

#endif