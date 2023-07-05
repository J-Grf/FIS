#include "multigrid.hpp"

size_t MG::gamma = 0;
size_t MG::nu1 = 0;
size_t MG::nu2 = 0;

MG::MG(const size_t _N, const size_t _l, const bool _wTO) : N{_N}, l{_l}, grid{getGrid(_N)}, writeToOutput{_wTO}, 
f_rhs{m_type(N + 1, std::vector<double>(N + 1, 0.0))} {
    std::cout << "N: " << N << " l: " << l << std::endl;
    for( size_t i = 1; i < N; i++) {
        for(size_t j = 1; j < N; j++) {
            f_rhs[i][j] = 8 * pow(M_PI, 2) * sin(2 * M_PI * grid[i * (N+1) + j].first) * sin(2 * M_PI * grid[i * (N+1) + j].second);
        }
    }

    //printFToFile();
}

void MG::printFToFile() {
    std::ofstream out;
    out.open("f.txt");
    for( size_t i = 0; i < N+1; i++) {
        for(size_t j = 0; j < N+1; j++) {
            out << f_rhs[i][j] << std::endl;
        }
    }
    out.close();
}

//maybe move
m_type MG::Restriction(const m_type& u, const size_t N_c) {
    std::ofstream out;
    if(writeToOutput)
        out.open("Restriction.txt");
    
    const size_t vecSize = N_c + 1;

    m_type f_coarse(vecSize, std::vector<double>(vecSize, 0.0));
    size_t ii{0}, jj{0};
    for(size_t i = 1; i < N_c; i++) {
        ii = 2*i;
        for(size_t j = 1; j < N_c; j++) {
            jj = 2*j;
            f_coarse[i][j] = N1D16 * (u[ii-1][jj-1] + 2.0 * u[ii][jj-1] + u[ii+1][jj-1] +
                                  2.0 * u[ii-1][jj] + 4.0 * u[ii][jj] + 2.0 * u[ii+1][jj] +
                                      u[ii-1][jj+1] + 2.0 * u[ii][jj+1] + u[ii+1][jj+1]);
        }
    }
    
    if(writeToOutput){
        for(size_t i = 0; i <= N_c; i++) {
            for(size_t j = 0; j <= N_c; j++) {
                out << f_coarse[i][j] << std::endl;
            }
        }
        out.close();
    }

    return f_coarse;
}

m_type MG::Prolongation(const m_type& u, const size_t N_c) {
    std::ofstream out;
    if(writeToOutput)
        out.open("Prolongation.txt");

    const size_t N = 2 * N_c;
    const size_t vecSize = N + 1;
    m_type f_fine;
    size_t ii{0}, jj{0};
    if(u.size() == 1) {
        m_type a(3, std::vector<double>(3, 0.0));
        f_fine = std::move(a);
        f_fine[1][1] = u[0][0];
    } else {
        m_type a(vecSize, std::vector<double>(vecSize, 0.0));
        f_fine = std::move(a);
        for(size_t i = 1; i < N_c; i++) {
            ii = 2*i;
            for(size_t j = 1; j < N_c; j++) {
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
    }

    if(writeToOutput) {
        for(size_t i = 0; i <= N; i++) {
            for(size_t j = 0; j <= N; j++) {
                out << f_fine[i][j] << std::endl;
            }
        }
        out.close();
    }

    return f_fine;
}
m_type MG::ComputeResidual(const m_type& f, const m_type& u, const size_t N) {
    const size_t vecSize = N + 1;
    m_type res(vecSize, std::vector<double>(vecSize, 0.0));
    const double N1DHSQ = pow(N, 2);

    for(size_t i = 1; i < N; i++) {
        for(size_t j = 1; j < N; j++) {
            res[i][j] = f[i][j] + ((u[i-1][j] - 2 * u[i][j] + u[i+1][j]) + (u[i][j-1] - 2 * u[i][j] + u[i][j+1])) * N1DHSQ;
        }
    }
    return res;
}

//recursive multigrid function
void MG::MG_Algorithm(const size_t l, m_type& u, const m_type& f) {
    std::cout << "MG_Algorithm: l=" << l << " gamma= " << gamma << " nu1= " << nu1 << " nu2= " << nu2 << std::endl;
    const size_t N = pow(2, l);
    const m_type u_tmp = GaussSeidel(u, f, nu1, N);
    const m_type res = ComputeResidual(f, u_tmp, N);

    m_type res_c = Restriction(res, N/2);

    m_type e;
    if(l == 1) { 
        std::cout << "reached l=1, N=" << N << std::endl; 
        //m_type zero(N+1, std::vector<double>(N+1, 0.0));
        m_type a(3, std::vector<double>(3, 0.0));
        e = std::move(a);
        e = GaussSeidel(e, -res, 1, N);
        e = {{e[1][1]}};
    } else {

        const size_t vec_c = pow(2, l-1) + 1;
        e.resize(vec_c);
        for(size_t i = 0; i < vec_c; i++) {
            e[i].resize(vec_c);
        }

        for(size_t j = 1; j <= gamma; j++) {
            //modifies e_coarse
            MG_Algorithm(l-1, e, -res_c);
        }
    }
  
    e = Prolongation(e, N/2);
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