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
    Matrix(const string _inputDir);
    size_t getDim() const { return dim; };

    void print() const;

friend void readFromFile (Matrix&);
friend vector<double> vectorProduct (const Matrix&, const std::vector<double>);
};

void readFromFile(Matrix& A);

std::vector<double> vectorProductCSR (const size_t dim, const std::vector<double>& IA, const std::vector<double>& J,
const std::vector<double>& V, const std::vector<double> x);

std::vector<double> vectorProductCSC (const size_t dim, const std::vector<double>& IA, const std::vector<double>& J,
const std::vector<double>& V, const std::vector<double> x);

vector<double> vectorProduct(const Matrix&, const std::vector<double>);

double dot(const std::vector<double>, const std::vector<double> );
