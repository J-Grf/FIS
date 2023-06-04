#ifndef _MATRIX_INCLUDE_
#define _MATRIX_INCLUDE_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

template<typename T>
const T eps = std::numeric_limits<T>::epsilon();

constexpr double Eps = 10E-8;
class Matrix {
    char sym_flag = -1;
    size_t dim = 0;
    size_t array_size = 0;
    
    std::vector<int> JM = {};
    std::vector<double> VM = {};

    std::string inputDir = "";

public:
    Matrix() = default;
    Matrix(const std::string& _inputDir);
    size_t getDim() const { return dim; };
    size_t getArrSize() const { return array_size; };
    double a_VM(const size_t i) const { return VM[i]; };
    int a_JM(const size_t i) const {return JM[i]; };
    std::vector<int> getJMVec() const { return JM; };
    std::vector<double> getVMVec() const { return VM; };

    void print() const;
    Matrix& operator=(const Matrix&) = default;

    friend void readFromFile (Matrix&);
    friend std::vector<double> vectorProduct (const Matrix&, const std::vector<double>&);
    friend std::vector<double> backwardSubMSR(const Matrix&, const std::vector<double>&, const size_t);
    friend std::vector<double> forwardSubMSR(const Matrix&, const std::vector<double>&, const size_t);
};

void readFromFile(Matrix& A);

//vectorProducts in different formats
std::vector<double> vectorProductCSR (const size_t dim, const std::vector<size_t>& IA, const std::vector<size_t>& J,
const std::vector<double>& V, const std::vector<double>& x);

std::vector<double> vectorProductCSC (const size_t dim, const std::vector<size_t>& IA, const std::vector<size_t>& J,
const std::vector<double>& V, const std::vector<double>& x);

std::vector<double> vectorProduct(const Matrix& A, const std::vector<double>& x);

std::vector<double> backwardSubMSR(const Matrix& A, const std::vector<double>& b, const size_t m);
std::vector<double> forwardSubMSR(const Matrix& A, const std::vector<double>& b, const size_t m);

//helper functions
inline double dotP(const std::vector<double> a, const std::vector<double> b){
    assert(b.size()==a.size());
    double res = 0;
    for (size_t i = 0; i < a.size(); i++) {
        res += a[i]*b[i];
    }
    return res;
}

inline double norm2(const std::vector<double> a) {
    return sqrt(dotP(a,a));
}

//conventional vectorProduct (but with transpose)
inline std::vector<double> VP(const std::vector<std::vector<double>>& A, const std::vector<double>& x, const size_t m_t) {
    assert(A[0].size() == x.size());

    std::vector<double> res(x.size());
    for(size_t i = 0; i < x.size(); i++) {
        for(size_t j = 0; j < m_t; j++) {
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
        v[i] -= a[i];
    return v;
}

template<typename T>
std::vector<T>& operator+=(std::vector<T>& v, const std::vector<T>& a) {
    assert(v.size() == a.size());
    for(size_t i = 0; i < v.size(); i++)
        v[i] += a[i];
    return v;
}

template<typename T>
std::vector<T> operator-(const std::vector<T>& v, const std::vector<T>& a) {
    assert(v.size() == a.size());
    std::vector<T> res(v.size());
    for(size_t i = 0; i < v.size(); i++)
        res[i] = v[i] - a[i];
    return res;
}

template<typename T>
std::vector<T> operator+(const std::vector<T>& v, const std::vector<T>& a) {
    assert(v.size() == a.size());
    std::vector<T> res(v.size());
    for(size_t i = 0; i < v.size(); i++)
        res[i] = v[i] + a[i];
    return res;
}

template<typename T>
std::vector<T> operator*(const T& a, const std::vector<T>& v) {
    std::vector<T> res(v.size());
    for(size_t i = 0; i < v.size(); i++) {
        res[i] = a * v[i];
    }
    return res;
}

#endif