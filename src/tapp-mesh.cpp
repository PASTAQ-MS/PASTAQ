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
    if (mData.size() == 0 || i > mDimensions.n - 1) {
        return std::nullopt;
    }

    auto deltaMz =
        (mBounds.maxMz - mBounds.minMz) / static_cast<double>(mDimensions.n - 1);
    return mBounds.minMz + deltaMz * i;
}

std::optional<double> Mesh::rt(unsigned int j) {
    if (mData.size() == 0 || j > mDimensions.m - 1) {
        return std::nullopt;
    }
    auto deltaRt =
        (mBounds.maxRt - mBounds.minRt) / static_cast<double>(mDimensions.m - 1);
    return mBounds.minRt + deltaRt * j;
}

// FIXME: Does this return a copy or a reference?
Grid::Dimensions Mesh::dim() { return mDimensions; }

Grid::Bounds Mesh::bounds() { return mBounds; }
