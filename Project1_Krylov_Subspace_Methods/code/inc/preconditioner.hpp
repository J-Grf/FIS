#ifndef _PRE_INCLUDE_
#define _PRE_INCLUDE_

#include "matrix.hpp"
#include <iostream>
#include <map>

enum PreConditioner {
    NONE,
    JACOBI,
    GAUSSSEIDEL,
    ILU
};

static std::map<std::string, PreConditioner> StringToPre { {"NONE", PreConditioner::NONE}, {"JACOBI", PreConditioner::JACOBI}, {"GAUSSSEIDEL", PreConditioner::GAUSSSEIDEL}, 
{"ILU", PreConditioner::ILU} };

inline std::ostream& operator<<(std::ostream& out, const PreConditioner value){
    const char* s = 0;
#define PROCESS_VAL(p) case(p): s = #p; break;
    switch(value){
        PROCESS_VAL(NONE);     
        PROCESS_VAL(JACOBI);     
        PROCESS_VAL(GAUSSSEIDEL);
        PROCESS_VAL(ILU);
    }
#undef PROCESS_VAL
    return out << s;
}

struct ILUout {   
    size_t nz, n;
    std::vector<double> PV;
    std::vector<int> PJM;
    std::vector<int> JD;

    ILUout(const size_t _nz, const size_t _n);
};

ILUout Ilu(const Matrix& A);
void applyPreConditioner(const Matrix& A, std::vector<double>& x, const PreConditioner PreCon);

#endif