#ifndef TAPP_GRID_MESH_HPP
#define TAPP_GRID_MESH_HPP

#include "tapp-grid.hpp"

class Mesh : public Grid::Interface {
   private:
    std::vector<double> mData;
    Grid::Dimensions mDimensions;
    Grid::Bounds mBounds;

   public:
    Mesh();
    ~Mesh();

    // Implementation methods for Grid::Interface.
    std::optional<double> at(int i, int j);
    std::optional<double> mz(int i);
    std::optional<double> rt(int j);
    Grid::Dimensions dim();
    Grid::Bounds bounds();
};

#endif /* TAPP_GRID_MESH_HPP */
