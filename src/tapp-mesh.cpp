#include <iostream>
#include <vector>

#include "tapp-mesh.hpp"

std::optional<double> Mesh::at(int i, int j) {
    return static_cast<double>(i) + static_cast<double>(j);
}
std::optional<double> Mesh::mz(int i) { return static_cast<double>(i); };
std::optional<double> Mesh::rt(int j) { return static_cast<double>(j); };
Grid::Dimensions Mesh::dim() { return mDimensions; }
Grid::Bounds Mesh::bounds() { return mBounds; }
