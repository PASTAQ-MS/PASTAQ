#include <cmath>
#include <iostream>

#include "tapp-weighted-mesh.hpp"

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
    mData[i + j * mDimensions.n] += value * weight;
    mWeights[i + j * mDimensions.n] += weight;
    ++mCounts[i + j * mDimensions.n];
    return true;
}

void WeightedMesh::printAll() {
    std::cout << "DATA:" << std::endl;
    for (unsigned int j = 0; j < mDimensions.m; ++j) {
        for (unsigned int i = 0; i < mDimensions.n; ++i) {
            auto e = at(i, j);
            if (e) {
                std::cout << e.value() << "\t";
            }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "WEIGHTS:" << std::endl;
    for (unsigned int j = 0; j < mDimensions.m; ++j) {
        for (unsigned int i = 0; i < mDimensions.n; ++i) {
            auto e = mWeights[i + j * mDimensions.n];
            std::cout << e << "\t";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
    std::cout << "COUNTS:" << std::endl;
    for (unsigned int j = 0; j < mDimensions.m; ++j) {
        for (unsigned int i = 0; i < mDimensions.n; ++i) {
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
    double deltaMz = (mBounds.maxMz - mBounds.minMz) /
                     static_cast<double>(mDimensions.n - 1);
    if (mz < mBounds.minMz || mz > mBounds.maxMz + deltaMz) {
        return std::nullopt;
    }
    double d = mz - mBounds.minMz;
    auto i = static_cast<unsigned int>(d / deltaMz);
    if (i > mDimensions.n - 1) {
        return std::nullopt;
    }
    return i;
}

std::optional<unsigned int> WeightedMesh::yIndex(double rt) {
    // In order to be consistent, the maximum value is rt + deltaRt. This
    // ensures all intervals contain the same number of points.
    double deltaRt = (mBounds.maxRt - mBounds.minRt) /
                     static_cast<double>(mDimensions.m - 1);
    if (rt < mBounds.minRt || rt > mBounds.maxRt + deltaRt) {
        return std::nullopt;
    }
    double d = rt - mBounds.minRt;
    auto j = static_cast<unsigned int>(d / deltaRt);
    if (j > mDimensions.m - 1) {
        return std::nullopt;
    }
    return j;
}

bool WeightedMesh::splash(double mzValue, double rtValue, double value) {
    // TODO: For now let's hardcode the magnitude of the splatting, will
    // parametrize later.
    double sigmaRt = 5.0;
    double sigmaMz = 100.0;

    // Get the gaussian square dimensions.
    double minRt = rtValue - 2 * sigmaRt;
    double maxRt = rtValue + 2 * sigmaRt;
    double minMz = mzValue - 2 * sigmaMz;
    double maxMz = mzValue + 2 * sigmaMz;

    // Even if the point lays outside the current grid, we still want to account
    // for it's contribution to the points in the frontier. However if the
    // minimum value in the splatting lays outside the grid, it will not have
    // any effect.
    auto iMin = xIndex(minMz);
    auto jMin = yIndex(minRt);
    if (iMin == std::nullopt && jMin == std::nullopt) {
        return false;
    }

    auto iMax = xIndex(maxMz) ? xIndex(maxMz) : mDimensions.n - 1;
    auto jMax = yIndex(maxRt) ? yIndex(maxRt) : mDimensions.m - 1;

    double x0 = mzValue;
    double y0 = rtValue;
    for (unsigned int j = jMin.value(); j <= jMax.value(); ++j) {
        for (unsigned int i = iMin.value(); i <= iMax.value(); ++i) {
            // No need to do boundary check, since we are sure we are inside the
            // grid.
            double x = mz(i).value();
            double y = rt(j).value();

            // Calculate the gaussian weight for this point.
            double a = (x - x0) / sigmaMz;
            double b = (y - y0) / sigmaRt;
            a *= a;
            b *= b;
            double weight = std::exp(-0.5 * (a + b));

            set(i, j, value * weight, weight);
            // TODO(alex): In this case the gaussian kernel is a square, se
            // could filter it to be a circle in exchange of extra computations.
            std::cout << "(" << i << "," << j << "," << weight << ")\t";
        }
        std::cout << std::endl;
    }
    return true;
}
