#include "tapp-mesh.hpp"

Mesh::Mesh(Grid::Dimensions dimensions, Grid::Bounds bounds)
    : mData(dimensions.n * dimensions.m),
      mDimensions(dimensions),
      mBounds(bounds) {}

std::optional<double> Mesh::at(unsigned int i, unsigned int j) {
    if (mData.size() == 0 || i > mDimensions.n - 1 || j > mDimensions.m - 1) {
        return std::nullopt;
    }
    return mData[i + j * mDimensions.m];
}

std::optional<double> Mesh::mz(unsigned int i) {
    return static_cast<double>(i);
}

std::optional<double> Mesh::rt(unsigned int j) {
    return static_cast<double>(j);
}

Grid::Dimensions Mesh::dim() { return mDimensions; }

Grid::Bounds Mesh::bounds() { return mBounds; }
