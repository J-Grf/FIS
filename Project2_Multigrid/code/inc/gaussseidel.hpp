#ifndef _GS_INCLUDE_
#define _GS_INCLUDE_

#include <vector>
using m_type = std::vector<std::vector<double>>;
constexpr double eps = 1E-10;

// GaussSeidel with output
m_type GaussSeidel(m_type& u, const m_type& u_ex, const m_type& f, const size_t nu, const size_t N);

//GaussSeidel without output
m_type GaussSeidel(m_type& u, const m_type& f, const size_t nu, const size_t N);

#endif