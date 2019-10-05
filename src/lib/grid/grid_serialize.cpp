#include "grid/grid_serialize.hpp"
#include "utils/serialization.hpp"

bool Grid::Serialize::read_grid(std::istream &stream, Grid *grid) {
    Serialization::read_uint64(stream, &grid->n);
    Serialization::read_uint64(stream, &grid->m);
    Serialization::read_uint64(stream, &grid->k);
    Serialization::read_uint64(stream, &grid->t);
    uint8_t instrument_type = 0;
    Serialization::read_uint8(stream, &instrument_type);
    Serialization::read_double(stream, &grid->reference_mz);
    Serialization::read_double(stream, &grid->fwhm_mz);
    Serialization::read_double(stream, &grid->fwhm_rt);
    Serialization::read_double(stream, &grid->min_mz);
    Serialization::read_double(stream, &grid->max_mz);
    Serialization::read_double(stream, &grid->min_rt);
    Serialization::read_double(stream, &grid->max_rt);
    grid->data = std::vector<double>(grid->n * grid->m);
    grid->bins_mz = std::vector<double>(grid->n);
    grid->bins_rt = std::vector<double>(grid->m);
    for (size_t i = 0; i < (grid->n * grid->m); ++i) {
        Serialization::read_double(stream, &grid->data[i]);
    }
    for (size_t i = 0; i < grid->n; ++i) {
        Serialization::read_double(stream, &grid->bins_mz[i]);
    }
    for (size_t i = 0; i < grid->m; ++i) {
        Serialization::read_double(stream, &grid->bins_rt[i]);
    }
    return stream.good();
}

bool Grid::Serialize::write_grid(std::ostream &stream, const Grid &grid) {
    Serialization::write_uint64(stream, grid.n);
    Serialization::write_uint64(stream, grid.m);
    Serialization::write_uint64(stream, grid.k);
    Serialization::write_uint64(stream, grid.t);
    Serialization::write_uint8(stream, grid.instrument_type);
    Serialization::write_double(stream, grid.reference_mz);
    Serialization::write_double(stream, grid.fwhm_mz);
    Serialization::write_double(stream, grid.fwhm_rt);
    Serialization::write_double(stream, grid.min_mz);
    Serialization::write_double(stream, grid.max_mz);
    Serialization::write_double(stream, grid.min_rt);
    Serialization::write_double(stream, grid.max_rt);
    for (size_t i = 0; i < (grid.n * grid.m); ++i) {
        Serialization::write_double(stream, grid.data[i]);
    }
    for (size_t i = 0; i < grid.n; ++i) {
        Serialization::write_double(stream, grid.bins_mz[i]);
    }
    for (size_t i = 0; i < grid.m; ++i) {
        Serialization::write_double(stream, grid.bins_rt[i]);
    }
    return stream.good();
}
