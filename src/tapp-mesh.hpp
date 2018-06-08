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
    std::optional<double> value_at(unsigned int i, unsigned int j) override;
    bool set_value(unsigned int i, unsigned int j, double value) override;
    std::optional<double> mz_at(unsigned int i) override;
    std::optional<double> rt_at(unsigned int j) override;
    std::optional<unsigned int> x_index(double mz) override;
    std::optional<unsigned int> y_index(double rt) override;
    Grid::Dimensions dim() override;
    Grid::Bounds bounds() override;
};

#endif /* TAPP_GRID_MESH_HPP */
