#ifndef GRID_INCLUDE
#define GRID_INCLUDE

#include <vector>
#include <utility>
#include <fstream>

inline std::vector<std::pair<double, double>> getGrid(const double N, const bool writeGrid = false, const std::string& File = "") {
    const size_t NSize = static_cast<size_t>(N) + 1;
    std::vector<std::pair<double, double>> grid(NSize * NSize);
    std::ofstream out;
    if(writeGrid)
        out.open(File);
    //coordinates of all grid points assuming (0,1) x (0,1)
    for(size_t i = 0; i < NSize; i++){
        for(size_t j = 0; j < NSize; j++) {
            grid[i * NSize + j] = std::make_pair( i/N, j/N);
        }
        if(writeGrid)
            out << grid[i * NSize].first << std::endl;
    }
    if(writeGrid)
        out.close();

    return grid;
}

inline m_type getExactSolution(const size_t N, const std::vector<std::pair<double,double>>& grid) {
    m_type u_ex(N+1, std::vector<double>(N+1, 0.0));

    for( size_t i = 1; i < N; i++) {
        for(size_t j = 1; j < N; j++) {
            u_ex[i][j] = sin(2 * M_PI * grid[i * (N+1) + j].first) * sin(2 * M_PI * grid[i * (N+1) + j].second);
        }
    }
    return u_ex;
}

inline void printSolution(const m_type& u, const std::string& File = "") {
    assert(u.size() == u[0].size());
    
    const size_t size = u.size();
    std::ofstream out;
    out.open(File);
    for(size_t i = 0; i < size; i++) {
        for(size_t j = 0; j < size; j++) {
            out << u[i][j] << std::endl;
        }
    }
    out.close();

}

#endif