#include "grid_serialize.hpp"
#include "utils/serialization.hpp"

bool Grid::Serialize::read_parameters(std::istream &stream,
                                      Grid::Parameters *parameters) {
    Serialization::read_uint64(stream, &parameters->dimensions.n);
    Serialization::read_uint64(stream, &parameters->dimensions.m);
    Serialization::read_double(stream, &parameters->bounds.min_rt);
    Serialization::read_double(stream, &parameters->bounds.max_rt);
    Serialization::read_double(stream, &parameters->bounds.min_mz);
    Serialization::read_double(stream, &parameters->bounds.max_mz);
    Serialization::read_double(stream, &parameters->smoothing_params.mz);
    Serialization::read_double(stream, &parameters->smoothing_params.sigma_mz);
    Serialization::read_double(stream, &parameters->smoothing_params.sigma_rt);
    Serialization::read_uint8(stream, &parameters->instrument_type);
    Serialization::read_uint8(stream, &parameters->flags);
    return stream.good();
}

bool Grid::Serialize::write_parameters(std::ostream &stream,
                                       const Grid::Parameters &parameters) {
    Serialization::write_uint64(stream, parameters.dimensions.n);
    Serialization::write_uint64(stream, parameters.dimensions.m);
    Serialization::write_double(stream, parameters.bounds.min_rt);
    Serialization::write_double(stream, parameters.bounds.max_rt);
    Serialization::write_double(stream, parameters.bounds.min_mz);
    Serialization::write_double(stream, parameters.bounds.max_mz);
    Serialization::write_double(stream, parameters.smoothing_params.mz);
    Serialization::write_double(stream, parameters.smoothing_params.sigma_mz);
    Serialization::write_double(stream, parameters.smoothing_params.sigma_rt);
    Serialization::write_uint8(stream, parameters.instrument_type);
    Serialization::write_uint8(stream, parameters.flags);
    return stream.good();
}

bool Grid::Serialize::read_point(std::istream &stream, Grid::Point *point) {
    Serialization::read_double(stream, &point->mz);
    Serialization::read_double(stream, &point->rt);
    Serialization::read_double(stream, &point->value);
    return stream.good();
}

bool Grid::Serialize::write_point(std::ostream &stream,
                                  const Grid::Point &point) {
    Serialization::write_double(stream, point.mz);
    Serialization::write_double(stream, point.rt);
    Serialization::write_double(stream, point.value);
    return stream.good();
}

bool Grid::Serialize::read_mesh(std::istream &stream, Grid::Mesh *mesh) {
    Serialization::read_uint64(stream, &mesh->n);
    Serialization::read_uint64(stream, &mesh->m);
    Serialization::read_uint64(stream, &mesh->k);
    Serialization::read_uint64(stream, &mesh->t);
    uint8_t instrument_type = 0;
    Serialization::read_uint8(stream, &instrument_type);
    Serialization::read_double(stream, &mesh->reference_mz);
    Serialization::read_double(stream, &mesh->fwhm_mz);
    Serialization::read_double(stream, &mesh->fwhm_rt);
    Serialization::read_double(stream, &mesh->min_mz);
    Serialization::read_double(stream, &mesh->max_mz);
    Serialization::read_double(stream, &mesh->min_rt);
    Serialization::read_double(stream, &mesh->max_rt);
    mesh->data = std::vector<double>(mesh->n * mesh->m);
    mesh->bins_mz = std::vector<double>(mesh->n);
    mesh->bins_rt = std::vector<double>(mesh->m);
    for (size_t i = 0; i < (mesh->n * mesh->m); ++i) {
        Serialization::read_double(stream, &mesh->data[i]);
    }
    for (size_t i = 0; i < mesh->n; ++i) {
        Serialization::read_double(stream, &mesh->bins_mz[i]);
    }
    for (size_t i = 0; i < mesh->m; ++i) {
        Serialization::read_double(stream, &mesh->bins_rt[i]);
    }
    return stream.good();
}

bool Grid::Serialize::write_mesh(std::ostream &stream, const Grid::Mesh &mesh) {
    Serialization::write_uint64(stream, mesh.n);
    Serialization::write_uint64(stream, mesh.m);
    Serialization::write_uint64(stream, mesh.k);
    Serialization::write_uint64(stream, mesh.t);
    Serialization::write_uint8(stream, mesh.instrument_type);
    Serialization::write_double(stream, mesh.reference_mz);
    Serialization::write_double(stream, mesh.fwhm_mz);
    Serialization::write_double(stream, mesh.fwhm_rt);
    Serialization::write_double(stream, mesh.min_mz);
    Serialization::write_double(stream, mesh.max_mz);
    Serialization::write_double(stream, mesh.min_rt);
    Serialization::write_double(stream, mesh.max_rt);
    for (size_t i = 0; i < (mesh.n * mesh.m); ++i) {
        Serialization::write_double(stream, mesh.data[i]);
    }
    for (size_t i = 0; i < mesh.n; ++i) {
        Serialization::write_double(stream, mesh.bins_mz[i]);
    }
    for (size_t i = 0; i < mesh.m; ++i) {
        Serialization::write_double(stream, mesh.bins_rt[i]);
    }
    return stream.good();
}
