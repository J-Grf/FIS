#ifndef GRID_INCLUDE
#define GRID_INCLUDE

#include <vector>
#include <utility>

inline std::vector<std::pair<double, double>> getGrid(const size_t N) {
    std::vector<std::pair<double, double>> grid(N);
    //coordinates of all grid points assuming (0,1) x (0,1)
    for(size_t i = 0; i < N; i++){
        for(size_t j = 0; j < N; j++) {
            grid.push_back(std::make_pair(i/N, j/N)); 
        }
    }
    return grid;
}

#endif