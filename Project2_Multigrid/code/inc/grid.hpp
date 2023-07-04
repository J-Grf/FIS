#ifndef GRID_INCLUDE
#define GRID_INCLUDE

#include <vector>
#include <utility>
#include <fstream>

inline std::vector<std::pair<double, double>> getGrid(const double N, const bool writeGrid = false, const std::string& File = "") {
    const size_t NSize = static_cast<size_t>(N);
    std::vector<std::pair<double, double>> grid(NSize * NSize);
    std::ofstream out;
    if(writeGrid)
        out.open(File);
    //coordinates of all grid points assuming (0,1) x (0,1)
    for(size_t i = 0; i < NSize; i++){
        for(size_t j = 0; j < NSize; j++) {
            grid[i * NSize + j] = std::make_pair((double)i/(N - 1.0), (double)j/(N - 1.0));
        }
        if(writeGrid)
            out << grid[i * NSize].first << std::endl;
    }
    if(writeGrid)
        out.close();

    return grid;
}

inline m_type getExactSolution(const size_t N, const std::vector<std::pair<double,double>>& grid) {
    m_type u_ex(N, std::vector<double>(N, 0.0));

    for( size_t i = 1; i < N - 1; i++) {
        for(size_t j = 1; j < N - 1; j++) {
            u_ex[i][j] = sin(2 * M_PI * grid[i * N + j].first) * sin(2 * M_PI * grid[i * N + j].second);
        }
    }
    return u_ex;
}

inline void printExactSolution(const m_type& u_ex, const std::string& File = "") {
    assert(u_ex.size() == u_ex[0].size());
    
    const size_t N = u_ex.size();
    std::ofstream out;
    out.open(File);
    for(size_t i = 0; i < N; i++) {
        for(size_t j = 0; j < N; j++) {
            out << u_ex[i][j] << std::endl;
        }
    }
    out.close();

}

#endif