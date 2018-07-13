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
