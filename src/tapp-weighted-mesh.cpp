#include "tapp-weighted-mesh.hpp"

#include <iostream>
WeightedMesh::WeightedMesh(Grid::Dimensions dimensions, Grid::Bounds bounds)
    : Mesh(dimensions, bounds),
      mWeights(dimensions.n * dimensions.m),
      mCounts(dimensions.n * dimensions.m) {}

bool WeightedMesh::set(unsigned int i, unsigned int j, double value,
                       double weight) {
    if (mData.size() == 0 || i > mDimensions.n - 1 || j > mDimensions.m - 1) {
        return false;
    }
    mData[i + j * mDimensions.m] += value * weight;
    mWeights[i + j * mDimensions.m] += weight;
    ++mCounts[i + j * mDimensions.m];
    return true;
}

void WeightedMesh::printAll() {
    std::cout << "DATA:" << std::endl;
    for (unsigned int i = 0; i < mDimensions.n; ++i) {
        for (unsigned int j = 0; j < mDimensions.m; ++j) {
            auto e = at(i, j);
            if (e) {
                std::cout << e.value() << "\t";
            }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "WEIGHTS:" << std::endl;
    for (unsigned int i = 0; i < mDimensions.n; ++i) {
        for (unsigned int j = 0; j < mDimensions.m; ++j) {
            auto e = mWeights[i + j * mDimensions.n];
            std::cout << e << "\t";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
    std::cout << "COUNTS:" << std::endl;
    for (unsigned int i = 0; i < mDimensions.n; ++i) {
        for (unsigned int j = 0; j < mDimensions.m; ++j) {
            auto e = mCounts[i + j * mDimensions.n];
            std::cout << e << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
