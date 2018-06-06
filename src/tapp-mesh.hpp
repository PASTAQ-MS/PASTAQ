#ifndef TAPP_GRID_MESH_HPP
#define TAPP_GRID_MESH_HPP

#include <vector>

#include "tapp-grid.hpp"

// This class represents regular mesh, mapped at discrete intervals for both mz
// and rt.
class Mesh : public Grid::Interface {
   protected:
    std::vector<double> mData;
    Grid::Dimensions mDimensions;
    Grid::Bounds mBounds;

   public:
    Mesh(Grid::Dimensions dimensions = {}, Grid::Bounds bounds = {});

    // Implementation methods for Grid::Interface.
    std::optional<double> at(unsigned int i, unsigned int j);
    std::optional<double> mz(unsigned int i);
    std::optional<double> rt(unsigned int j);
    Grid::Dimensions dim();
    Grid::Bounds bounds();
};

#endif /* TAPP_GRID_MESH_HPP */
