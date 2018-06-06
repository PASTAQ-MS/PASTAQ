#include "tapp-weighted-mesh.hpp"

#include <iostream>
WeightedMesh::WeightedMesh(Grid::Dimensions dimensions, Grid::Bounds bounds)
    : Mesh(dimensions, bounds),
      mWeights(dimensions.n * dimensions.m),
      mCounts(dimensions.n * dimensions.m) {}

std::optional<double> WeightedMesh::weightAt(unsigned int i, unsigned int j) {
    if (mData.size() == 0 || i > mDimensions.n - 1 || j > mDimensions.m - 1) {
        return std::nullopt;
    }
    return mWeights[i + j * mDimensions.m];
}

std::optional<double> WeightedMesh::countsAt(unsigned int i, unsigned int j) {
    if (mData.size() == 0 || i > mDimensions.n - 1 || j > mDimensions.m - 1) {
        return std::nullopt;
    }
    return mCounts[i + j * mDimensions.m];
}

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

std::optional<unsigned int> WeightedMesh::xIndex(double mz) {
    // In order to be consistent, the maximum value is mz + deltaMz. This
    // ensures all intervals contain the same number of points.
    double deltaMz = (mBounds.maxMz - mBounds.minMz) / mDimensions.n;
    if (mz < mBounds.minMz || mz > mBounds.maxMz + deltaMz) {
        return std::nullopt;
    }
    double d = mz - mBounds.minMz;
    auto i = static_cast<unsigned int>(d / deltaMz);
    if (i > mDimensions.n) {
        return std::nullopt;
    }
    return i;
}

std::optional<unsigned int> WeightedMesh::yIndex(double rt) {
    // In order to be consistent, the maximum value is rt + deltaRt. This
    // ensures all intervals contain the same number of points.
    double deltaRt = (mBounds.maxRt - mBounds.minRt) / mDimensions.m;
    if (rt < mBounds.minRt || rt > mBounds.maxRt + deltaRt) {
        return std::nullopt;
    }
    double d = rt - mBounds.minRt;
    auto j = static_cast<unsigned int>(d / deltaRt);
    if (j > mDimensions.m) {
        return std::nullopt;
    }
    return j;
}
