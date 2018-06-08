#ifndef TAPP_GRID_MESH_HPP
#define TAPP_GRID_MESH_HPP

#include <vector>

#include "tapp-grid.hpp"

// This class represents regular mesh, mapped at discrete intervals for both mz
// and rt.
class Mesh : public Grid::Interface {
   protected:
    std::vector<double> m_data;
    Grid::Dimensions m_dimensions;
    Grid::Bounds m_bounds;

   public:
    Mesh(Grid::Dimensions dimensions = {}, Grid::Bounds bounds = {});

    // Implementation methods for Grid::Interface.
    std::optional<double> at(unsigned int i, unsigned int j);
    bool set(unsigned int i, unsigned int j, double value);
    std::optional<double> mz(unsigned int i);
    std::optional<double> rt(unsigned int j);
    Grid::Dimensions dim();
    Grid::Bounds bounds();
};

#endif /* TAPP_GRID_MESH_HPP */
