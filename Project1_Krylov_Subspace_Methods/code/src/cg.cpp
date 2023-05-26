
#include <utility>
#include "matrix.hpp"
#include "cg.hpp"

// Conjugate Gradient method only for s.p.d matrices!
std::pair<std::vector<double>, std::vector<double>> CGMethod(const Matrix& A, const std::vector<double>& b, 
const std::vector<double>& x0, const size_t m) {
    const std::vector<double> r0 = b - vectorProduct(A, x0);
    std::vector<double> r = r0, p = r, x = x0, r_old = r;

    double alpha = 0.0, beta = 0.0, relRes = 0.0;
    for(size_t i = 0; i < m; i++) {
        alpha = dotP(r, r) / dotP(vectorProduct(A,p), p);
        x += alpha * p;
        r_old = r;
        r -= alpha * vectorProduct(A,p);    
        beta = dotP(r, r) / dotP(r_old, r_old);
        p = r + beta * p;


        //TODO: Maybe comnpute spectral radius?
        relRes = norm2(r) / norm2(r0);

        printCG(x, r, p, relRes, i);
    }

    return std::make_pair(x, r);

}

void printCG(const std::vector<double>& x, const std::vector<double>& r, const std::vector<double>& p, double relRes, const size_t it) {
    assert(x.size() == r.size() && x.size() == p.size());

    std::string OUTS;
    std::cout << "CG Iteration " << it << std::endl;
    OUTS += "x : [";
    for(size_t i = 0; i < x.size(); i++) {
        OUTS += std::to_string(x[i]) + ", ";
    }
    OUTS += "] \n";

    OUTS += "r : [";
    for(size_t i = 0; i < x.size(); i++) {
        OUTS += std::to_string(r[i]) + ", ";
    }
    OUTS += "] ";
    OUTS += "rel res: " + std::to_string(relRes);
    OUTS += "\n";

    OUTS += "p : [";
    for(size_t i = 0; i < x.size(); i++) {
        OUTS += std::to_string(p[i]) + ", ";
    }
    OUTS += "] \n";

    std::cout << OUTS << std::endl;
}
