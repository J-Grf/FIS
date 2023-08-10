#ifndef _EVMETHODS_INCLUDE_
#define _EVMETHODS_INCLUDE_

#include <vector>
#include <memory>
#include "matrix.hpp"

constexpr double EPS = 1E-12;
constexpr double lambdaCG = 9.5986080894852857E3;

struct PowItObj{
    size_t maxIdx{};
    const int bufSize = 1E5;
    double lambdaMax{};
    std::vector<double> diff;
    std::vector<double> times;
    std::vector<double> error;

    PowItObj();
};

double powerIteration(const Matrix& A, const std::vector<double>& q0, std::unique_ptr<PowItObj>& PowItPtr, double eps = 1E-8);

double LanczosMethod(const Matrix& A, std::vector<double>& v, const size_t m, double customEps = -1.0);

void printAlphaBeta(const std::vector<double>& alpha, const std::vector<double>& beta);

#endif