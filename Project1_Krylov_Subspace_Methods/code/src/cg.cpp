
#include <utility>
#include "matrix.hpp"
#include "cg.hpp"

// Conjugate Gradient method only for s.p.d matrices!
std::pair<std::vector<double>, std::vector<double>> CGMethod(const Matrix& A, const std::vector<double>& b, 
const std::vector<double>& x0, const size_t m) {
    const std::vector<double> r0 = b - vectorProduct(A, x0);
    std::vector<double> r = r0, p = r, x = x0, r_old = r, error(x0.size()), x_ex(x0.size(),1.0), eANorm, rNorm;

    double alpha = 0.0, beta = 0.0, relRes = 0.0;
    for(size_t i = 0; i < m; i++) {
        std::cout << "CG Iteration " << i + 1 << std::endl;
        alpha = dotP(r, r) / dotP(vectorProduct(A,p), p);
        x += alpha * p;
        error = x - x_ex;
        eANorm.push_back(sqrt(dotP(vectorProduct(A,error), error)));

        r_old = r;
        r -= alpha * vectorProduct(A,p);    
        beta = dotP(r, r) / dotP(r_old, r_old);
        p = r + beta * p;

        rNorm.push_back(norm2(r));
        relRes = norm2(r) / norm2(r0);
        if(relRes < Eps)
            break;

        //printCG(x, r, p, relRes, i);
    }
    saveData(eANorm, rNorm);

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

void saveData(const std::vector<double>& eANorm, const std::vector<double>& rNorm) {

    std::ofstream out;
    out.open("eANorm.txt");
    for(size_t i = 0; i < eANorm.size(); i++) {
        out << eANorm[i] << std::endl;
    }
    out.close();

    out.open("rNorm.txt");
    for(size_t i = 0; i < rNorm.size(); i++) {
        out << rNorm[i] << std::endl;
    }
    out.close();
}
