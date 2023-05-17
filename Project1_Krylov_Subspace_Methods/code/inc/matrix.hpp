#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

class Matrix {
    char sym_flag = -1;
    size_t dim = 0;
    size_t array_size = 0;
    
    vector<int> JM = {};
    vector<double> VM = {};

    const string inputDir;

public:
    Matrix(const string);

    void print() const;

friend void readFromFile (Matrix&);
friend vector<double> vectorProductCSR (const Matrix&, const std::vector<double>);
friend vector<double> vectorProductCSC (const Matrix&, const std::vector<double>);
friend vector<double> vectorProduct (const Matrix&, const std::vector<double>);
};

void readFromFile(Matrix&);

vector<double> vectorProductCSR(const Matrix&, const std::vector<double>);
vector<double> vectorProductCSC(const Matrix&, const std::vector<double>);
vector<double> vectorProduct(const Matrix&, const std::vector<double>);