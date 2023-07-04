#include "multigrid.hpp"

size_t MG::gamma = 0;
size_t MG::nu1 = 0;
size_t MG::nu2 = 0;

MG::MG(const size_t _N, const size_t _l) : N{_N}, l{_l}, f{m_type(N, std::vector<double>(N, 0.0))}, grid{getGrid(_N)} {
    std::cout << "N: " << N << " l: " << l << std::endl;
    for( size_t i = 1; i < N - 1; i++) {
        for(size_t j = 1; j < N - 1; j++) {
            f[i][j] = 8 * pow(M_PI, 2) * sin(2 * M_PI * grid[i * N + j].first) * sin(2 * M_PI * grid[i * N + j].second);
        }
    }
}

//maybe move
m_type MG::Restriction(const m_type& u, const size_t N_c) {
    std::ofstream out;
    out.open("Restriction.txt");

    m_type f_coarse(N_c, std::vector<double>(N_c, 0.0));
    size_t ii{0}, jj{0};
    for(size_t i = 1; i < N_c - 1; i++) {
        ii = 2*i;
        for(size_t j = 1; j < N_c - 1; j++) {
            jj = 2*j;
            f_coarse[i][j] = N1D16 * (u[ii-1][jj-1] + 2.0 * u[ii][jj-1] + u[ii+1][jj-1] +
                                  2.0 * u[ii-1][jj] + 4.0 * u[ii][jj] + 2.0 * u[ii+1][jj] +
                                      u[ii-1][jj+1] + 2.0 * u[ii][jj+1] + u[ii+1][jj+1]);
        }
    }
    
    for(size_t i = 0; i < N_c; i++) {
        for(size_t j = 0; j < N_c; j++) {
            out << f_coarse[i][j] << std::endl;
        }
    }
    out.close();
    return f_coarse;
}

m_type MG::Prolongation(const m_type& u, const size_t N_c) {
    std::ofstream out;
    out.open("Prolongation.txt");

    const size_t N = 2 * N_c;
    m_type f_fine(N, std::vector<double>(N, 0.0));
    size_t ii{0}, jj{0};
    for(size_t i = 1; i < N_c - 1; i++) {
        ii = 2*i;
        for(size_t j = 1; j < N_c - 1; j++) {
            jj = 2*j;
            f_fine[ii-1][jj-1] += N1D4 * u[i][j];
            f_fine[ii][jj-1] += N1D2 * u[i][j];
            f_fine[ii+1][jj-1] += N1D4 * u[i][j];
            
            f_fine[ii-1][jj] += N1D2 * u[i][j];
            f_fine[ii][jj] += u[i][j];
            f_fine[ii+1][jj] += N1D2 * u[i][j];
            
            f_fine[ii-1][jj+1] += N1D4 * u[i][j];
            f_fine[ii][jj+1] += N1D2 * u[i][j];
            f_fine[ii+1][jj+1] += N1D4 * u[i][j];
        }
    }
    
    for(size_t i = 0; i < N; i++) {
        for(size_t j = 0; j < N; j++) {
            out << f_fine[i][j] << std::endl;
        }
    }
    out.close();

    return f_fine;
}
m_type MG::ComputeResidual(const m_type& f, const m_type& u) {
    m_type res(N, std::vector<double>(N, 0.0));
    const double N1DHSQ = pow((N-1), 2);

    for(size_t i = 1; i < N - 1; i++) {
        for(size_t j = 1; j < N - 1; j++) {
            res[i][j] = f[i][j] + (u[i-1][j] - 2 * u[i][j] + u[i+1][j]) * N1DHSQ + (u[i][j-1] - 2 * u[i][j] + u[i][j+1]) * N1DHSQ;
        }
    }
    return res;
}

//recursive multigrid function
void MG::MG_Algorithm(const size_t l, m_type& u, const m_type& f) {
    const m_type u_tmp = GaussSeidel(u, f, nu1, N);
    const m_type res = ComputeResidual(f, u_tmp);

    m_type res_c = Restriction(res, N/2);

    m_type e_coarse;
    if(l == 1) {
        const m_type res_neg{{-res_c[0][0]}}; 
        m_type zero{{0.0}};
        m_type e_0 = GaussSeidel(zero, res_neg, 1, 2);
    } else {
        for(size_t j = 1; j <= gamma; j++) {
            //implement -operator
            MG_Algorithm(l-1, e_coarse, -res_c);
        }
    }

    m_type e = Prolongation(e_coarse, N);
    //implement operator - for m_type
    m_type diff = u_tmp - e;

    u = GaussSeidel(diff, f, nu2, N);
}

//unary minus
m_type operator-(const m_type& v) {
    m_type u = v;
    for(size_t i = 0; i < v.size(); i++) {
        for(size_t j = 0; j < v[0].size(); j++){
            u[i][j] = -v[i][j];
        }
    }
    return u;
}

//binary minus
m_type operator-(const m_type& u, const m_type& v) {
    assert(u.size() == v.size());
    m_type r = u;
    for(size_t i = 0; i < v.size(); i++) {
        for(size_t j = 0; j < v[0].size(); j++){
            r[i][j] = u[i][j] - v[i][j];
        }
    }
    return r;
}