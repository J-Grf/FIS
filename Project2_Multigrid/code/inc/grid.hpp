#ifndef GRID_INCLUDE
#define GRID_INCLUDE

#include <vector>
#include <utility>
#include <fstream>

inline std::vector<std::pair<double, double>> getGrid(const double N, const bool writeGrid = false) {
    const size_t NSize = static_cast<size_t>(N);
    std::vector<std::pair<double, double>> grid(NSize * NSize);
    std::ofstream out;
    if(writeGrid)
        out.open("grid.txt");
    //coordinates of all grid points assuming (0,1) x (0,1)
    for(size_t i = 0; i < NSize; i++){
        for(size_t j = 0; j < NSize; j++) {
            grid[i * NSize + j] = std::make_pair(i/(N - 1), j/(N - 1));
        }
        if(writeGrid)
            out << grid[i * NSize].first << std::endl;
    }
    if(writeGrid)
        out.close();

    return grid;
}

#endif