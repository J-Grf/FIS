#include "gmres.hpp"
#include <boost/timer/timer.hpp>

// regular backwardSubstitution 
std::vector<double> backwardSub(const matrixType<double>& A, const std::vector<double>& b, const int m) {
    //use strict upper matrix of Hessenberg matrix
    assert(A[0].size() == b.size() - 1 && "Dimensions of A and b do not coincide, backwardSub not possible!");
    
    std::vector<double> x(m, 0);
    int idx = m - 1;

    x[idx] = b[idx] / A[idx][idx];
    for(int i = idx - 1 ; i >= 0; i--) {
        x[i] = b[i];
        for(int j = i+1; j < m; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }

    return x;
}

//Minimal Residual method
std::vector<double> MR_method(const Matrix& A, const std::vector<double>& b, const std::vector<double>& x0) {
    using namespace std;

    vector<double> r0 = b - vectorProduct(A, x0);
    vector<double> p,r,x(A.getDim());
    r = r0;
    double alpha;

    while( norm2(r)/norm2(r0) > Eps) {

        p = vectorProduct(A, r);
        alpha = dotP(r, p) / dotP(p, p);
        x += alpha * r;
        r -= alpha * p;
        std::cout << "relR for MR " << norm2(r)/norm2(r0) << std::endl;
        break;
    } 

    return x;
}

// GMRES
std::pair<std::vector<double>, double> GMRES(const Matrix& A, const std::vector<double>& x0, const std::vector<double>& b, const size_t m, 
const PreConditioner PreCon = NONE) {
    const size_t num = m + 1;

    std::vector<double> r0 = b - vectorProduct(A, x0);
    matrixType<double> V(1);
    std::unique_ptr<ILUout> ILUobj;

#ifndef DISABLEIO
    if(PreCon != NONE) {
        // left-preconditioning: in this case, r0 will be rbar -> v1
        applyPreConditioner(A, r0, PreCon, ILUobj);
    }
#endif

    const double normR0 = norm2(r0);
    const double invBeta = 1.0 / normR0; 
    
    // v1
    for(size_t i = 0; i < r0.size(); i++) {
        V[0].push_back(invBeta * r0[i]);
    }
    
    // H_u: upper Hessenberg 
    matrixType<double> H;
    for(size_t i = 0; i < num; i++) {
        H.push_back(std::vector<double>(m, 0));
    }

    std::vector<double> g(num, 0.0);
    g[0] = normR0;

#ifndef DISABLEIO
    std::cout << "norm2(r0): " << g[0] << std::endl; 
#endif

    std::vector<double> c, s, hj, relRes;
    relRes.push_back(1.0);
    size_t j;
    for(j = 0; j < m ; j++) {
    
#ifndef DISABLEIO
        std::cout << "-------GMRES sub-iteration " << j << " --------" <<  std::endl;
#endif
        hj = getKrylov(A, V, H, j, PreCon, ILUobj);
    
        /* for(size_t i = 0; i < hj.size(); i++) {
            std::cout << "hj[" << i << "]: " << hj[i] << std::endl;
        }
        
        //print upper Hessenberg
        for(size_t i = 0; i < H.size(); i++) {
            for(size_t j = 0; j < H[0].size(); j++){
                std::cout << "H[" << i << "][" << j << "]: " << H[i][j] << std::endl;
            }
        }
        printKrylov(V); */

        //Apply previous rotations
        double tmp;
        for(size_t k = 1; k <= j; k++) {
            tmp = hj[k - 1]; // store before overwriting
            hj[k - 1] = c[k-1] * tmp + s[k-1] * hj[k];
            H[k-1][j] = hj[k - 1];
            
            hj[k] = -s[k-1] * tmp + c[k-1] * hj[k];
            H[k][j] = hj[k];
        }

        //new rotations
        const double alpha = sqrt(hj[j] * hj[j] + hj[j+1] * hj[j+1]);
        c.push_back(hj[j] / alpha);
        s.push_back(hj[j+1] / alpha);
        hj[j] = alpha;
        H[j][j] = alpha;

        g[j+1] = -s[j] * g[j];
        g[j] *= c[j];
        //premature exit
        relRes.push_back(std::abs(g[j+1] / normR0));
        if(relRes.back() < Eps) {
            std::cout << "exiting with rel residual "  << relRes.back() << std::endl;
            break;
        }
    }
    
    //write rel residuals to file
#ifndef DISABLEIO
    saveRelResiduals(relRes, PreCon);
    if(PreCon == NONE)
        saveDotPofKrylovVectors(V);
#endif
    
    size_t m_tilde = std::min(j + 1 , m);
#ifndef DISABLEIO
    std::cout << "m_tilde: " << m_tilde << std::endl;
#endif
    std::vector<double> xm, y;

    y = backwardSub(H, g, m_tilde);
    for(size_t i = m_tilde; i < x0.size(); i++) {
        y.push_back(0);
    }

    xm = x0;
    std::vector<double> tmp = VP(V,y, m_tilde);
    for(size_t i = 0; i < xm.size(); i++) {
        xm[i] += tmp[i];
    }

    double rnorm = std::abs(g[m_tilde]);
    return make_pair(xm, rnorm);
}

// restarted GMRES
std::vector<double> GMRES_Res(const Matrix& A, const std::vector<double>& x0, const std::vector<double>& b, const size_t m, const PreConditioner PreCon = NONE){
    assert(A.getDim() == b.size());
    using namespace std;
    vector<double> r0 = b - vectorProduct(A, x0);
    const double r0Norm  = norm2(r0);

    size_t it = 0;
    double rho = 1.0;
    vector<double> x = x0;
    vector<double> r;
    const size_t restartPar = m;

    boost::timer::cpu_timer timer;

    std::cout << "Restarted GMRES with Preconditioner " << PreCon << " and " << m << " Krylov Vectors" << std::endl;
    
    timer.start();
    while(rho > Eps) {
        it++;
        if(it == 2)
            std::cout << "restart" << std::endl;

#ifndef DISABLEIO
        std::cout << "-----------------Iteration number " << it << " -------------------" << std::endl;
#endif
        const pair<vector<double>, double> res = GMRES(A, x, b, restartPar, PreCon);
        x = res.first;
        rho = res.second / r0Norm;

        // save residual
#ifndef DISABLEIO
        r.push_back(res.second);
        std::cout << "residual " << r.back() << std::endl;
        std::cout << "rel residual " << rho << std::endl;
        /*for(size_t i = 0; i < x.size(); i++){
            std::cout << "x[" << i << "]: " << x[i] << std::endl;
        }*/
#endif

        // for comparison with MR and GMRES(1)
        //if(m == 1)
        //    break;
        
    }
    timer.stop();
    std::cout  << "Restarted GMRES thread timer: " << timer.format() << std::endl
               << "m = " << m << " iterations = " << it << std::endl;

    return x;
}

void saveRelResiduals(const std::vector<double>& relRes, const PreConditioner PreCon) {
    std::string PreConStr;
    for(const auto& p : StringToPre) {
        if(p.second == PreCon) {
            PreConStr = p.first;
        }
    }

    std::ofstream out;
    out.open("relResiduals_" + PreConStr + ".txt");
    for(size_t i = 0; i < relRes.size(); i++) {
        out << relRes[i] << std::endl;
    }
    out.close();
}

void saveDotPofKrylovVectors(const matrixType<double>& V) {
    std::ofstream out;
    out.open("DotPKrylov.txt");
    for(size_t i = 0; i < V.size(); i++) {
        out << dotP(V[0], V[i]) << std::endl;
    }
    out.close();
}