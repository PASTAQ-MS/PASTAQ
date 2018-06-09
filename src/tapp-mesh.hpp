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
    Instrument::Config m_instrument_config;

   public:
    Mesh(Grid::Dimensions dimensions = {}, Grid::Bounds bounds = {},
         Instrument::Config instrument_config = {Instrument::QUAD, 200, 0.001});

    // Implementation methods for Grid::Interface.
    std::optional<double> value_at(unsigned int i, unsigned int j) override;
    bool set_value(unsigned int i, unsigned int j, double value) override;
    std::optional<double> mz_at(unsigned int i) override;
    std::optional<double> rt_at(unsigned int j) override;
    std::optional<unsigned int> x_index(double mz) override;
    std::optional<unsigned int> y_index(double rt) override;
    double sigma_at_mz(double mz) override;
    Grid::Dimensions dim() override;
    Grid::Bounds bounds() override;
};

#endif /* TAPP_GRID_MESH_HPP */
