#ifndef _MG_INCLUDE_
#define _MG_INCLUDE_

#include <iostream>
#include "gaussseidel.hpp"
#include "grid.hpp"

static constexpr double N1D16 = 0.0625; 
static constexpr double N1D4  = 0.25;
static constexpr double N1D2  = 0.5;

class MG {
    // indexing of grid goes from 0...N, so all loops over inner points
    // should go from i = 1 to i < N
    size_t N{0};
    size_t l{0};
    static size_t gamma;
    static size_t nu1;
    static size_t nu2;

    std::vector<std::pair<double, double>> grid;

    bool writeToOutput = false;

    public:
    m_type f_rhs;
    static void setStaticVariables(const size_t _g, const size_t _nu1, const size_t _nu2){
        gamma = _g;
        nu1 = _nu1;
        nu2 = _nu2;
    }
    
    MG(const size_t _N, const size_t _l, const bool _wTO);

    m_type Restriction(const m_type& r_fine, const size_t N_c);
    m_type Prolongation(const m_type& u_coarse, const size_t N_c);
    m_type ComputeResidual(const m_type& f, const m_type& u, const size_t N);
    void MG_Algorithm(const size_t l, m_type& u, const m_type& f);

    void printFToFile();
};

//unary minus
m_type operator-(const m_type& v);

//binary minus
m_type operator-(const m_type& u, const m_type& v);


#endif
