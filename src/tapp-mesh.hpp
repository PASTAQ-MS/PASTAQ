#ifndef TAPP_GRID_MESH_HPP
#define TAPP_GRID_MESH_HPP

#include <vector>

#include "tapp-grid.hpp"

// This class represents regular mesh, mapped at discrete intervals for both mz
// and rt.
class RegularMesh : public Grid::Interface {
   protected:
    std::vector<double> m_data;

    // Required parameters.
    Grid::Dimensions m_dimensions;
    Grid::Bounds m_bounds;
    Instrument::Type m_instrument_type;
    Grid::SmoothingParams m_smoothing_params;

   public:
    RegularMesh(Grid::Dimensions dimensions, Grid::Bounds bounds,
                Instrument::Type instrument_type,
                Grid::SmoothingParams smoothing_params);

    // Load the data from a istream.
    // TODO(alex): It would be a good idea to put all the parameters directly as
    // binary data instead of reading from two different files, .mesh and .dat.
    // The header on the .dat file will contain all the relevant parameters
    // necessary to interpret the data.
    //     see:
    //     https://stackoverflow.com/questions/416436/what-to-put-in-a-binary-data-files-header
    bool load_dat(std::istream &stream, Grid::Dimensions dimensions,
                  Grid::Bounds bounds, Instrument::Type instrument_type,
                  Grid::SmoothingParams smoothing_params);

    // Writes the current data to the stream.
    bool write_dat(std::ostream &stream);

    // Implementation methods for Grid::Interface.
    std::optional<double> value_at(uint32_t i, uint32_t j) override;
    bool set_value(uint32_t i, uint32_t j, double value) override;
    std::optional<double> mz_at(uint32_t i) override;
    std::optional<double> rt_at(uint32_t j) override;
    std::optional<uint32_t> x_index(double mz) override;
    std::optional<uint32_t> y_index(double rt) override;
    double sigma_mz(double mz) override;
    double sigma_rt() override;
    Grid::Dimensions dim() override;
    Grid::Bounds bounds() override;
};

#endif /* TAPP_GRID_MESH_HPP */
