#include <string>

#include "grid_files.hpp"
#include "utils/serialization.hpp"

bool Grid::Files::Dat::read(std::istream &stream,
                            std::vector<double> *destination,
                            Grid::Parameters *parameters) {
    if (destination == nullptr || parameters == nullptr) {
        return false;
    }
    if (Grid::Files::Dat::read_parameters(stream, parameters)) {
        auto n_points = parameters->dimensions.n * parameters->dimensions.m;
        destination->resize(n_points);
        stream.seekg(0, std::ios::beg);
        stream.read(reinterpret_cast<char *>(&(*destination)[0]),
                    sizeof(double) * n_points);
    } else {
        return false;
    }
    return stream.good();
}

bool Grid::Files::Dat::write(std::ostream &stream,
                             const std::vector<double> &source,
                             const Grid::Parameters &parameters) {
    stream.write(reinterpret_cast<const char *>(&source[0]),
                 sizeof(double) * source.size());
    Grid::Files::Dat::write_parameters(stream, parameters);
    return stream.good();
}

bool Grid::Files::Dat::read_range(std::istream &stream,
                                  const Grid::Bounds &bounds,
                                  std::vector<double> *destination,
                                  Grid::Parameters *parameters) {
    if (destination == nullptr || parameters == nullptr) {
        return false;
    }
    Grid::Parameters file_parameters = {};
    if (Grid::Files::Dat::read_parameters(stream, &file_parameters)) {
        // Get the indexes for the sliced range.
        auto i_min = Grid::x_index(bounds.min_mz, file_parameters);
        auto i_max = Grid::x_index(bounds.max_mz, file_parameters);
        auto j_min = Grid::y_index(bounds.min_rt, file_parameters);
        auto j_max = Grid::y_index(bounds.max_rt, file_parameters);

        // Check that indexes min/max are not swapped.
        if (i_min > i_max || j_min > j_max) {
            return false;
        }

        // Check that indexes are inside the grid.
        if ((i_min > file_parameters.dimensions.n - 1) ||
            (j_min > file_parameters.dimensions.m - 1)) {
            return false;
        }

        // Get the precise bounds for sliced range.
        Grid::Bounds sliced_bounds = {};
        sliced_bounds.min_mz = Grid::mz_at(i_min, file_parameters);
        sliced_bounds.max_mz = Grid::mz_at(i_max, file_parameters);
        sliced_bounds.min_rt = Grid::rt_at(j_min, file_parameters);
        sliced_bounds.max_rt = Grid::rt_at(j_max, file_parameters);

        // Get dimensions for sliced range.
        Grid::Dimensions sliced_dimensions = {};
        sliced_dimensions.n = (i_max - i_min) + 1;
        sliced_dimensions.m = (j_max - j_min) + 1;

        // Set up parameters.
        parameters->dimensions = sliced_dimensions;
        parameters->bounds = sliced_bounds;
        parameters->smoothing_params = file_parameters.smoothing_params;
        parameters->flags = file_parameters.flags;
        parameters->instrument_type = file_parameters.instrument_type;

        // Allocate memory on destination array.
        auto n_points = parameters->dimensions.n * parameters->dimensions.m;
        destination->resize(n_points);

        // Fill the destination with the data from the sliced parameters.
        size_t read_offset = 0;
        for (size_t j = j_min; j <= j_max; ++j) {
            stream.seekg(
                (j * file_parameters.dimensions.n + i_min) * sizeof(double),
                std::ios::beg);
            stream.read(reinterpret_cast<char *>(&(*destination)[read_offset]),
                        sizeof(double) * parameters->dimensions.n);
            read_offset += parameters->dimensions.n;
            if (stream.bad()) {
                return false;
            }
        }
    } else {
        return false;
    }
    return stream.good();
}

bool Grid::Files::Dat::write_range(std::ostream &stream,
                                   const Grid::Bounds &bounds,
                                   const std::vector<double> &source,
                                   const Grid::Parameters &parameters) {
    auto i_min = Grid::x_index(bounds.min_mz, parameters);
    auto i_max = Grid::x_index(bounds.max_mz, parameters);
    auto j_min = Grid::y_index(bounds.min_rt, parameters);
    auto j_max = Grid::y_index(bounds.max_rt, parameters);

    // Check that indexes min/max are not swapped.
    if (i_min > i_max || j_min > j_max) {
        return false;
    }

    // Check that indexes are inside the grid.
    if ((i_min > parameters.dimensions.n - 1) ||
        (j_min > parameters.dimensions.m - 1)) {
        return false;
    }

    // Get the precise bounds for sliced range.
    Grid::Bounds sliced_bounds = {};
    sliced_bounds.min_mz = Grid::mz_at(i_min, parameters);
    sliced_bounds.max_mz = Grid::mz_at(i_max, parameters);
    sliced_bounds.min_rt = Grid::rt_at(j_min, parameters);
    sliced_bounds.max_rt = Grid::rt_at(j_max, parameters);

    // Get dimensions for sliced range.
    Grid::Dimensions sliced_dimensions = {};
    sliced_dimensions.n = (i_max - i_min) + 1;
    sliced_dimensions.m = (j_max - j_min) + 1;

    // Set up parameters.
    Grid::Parameters sliced_parameters = {};
    sliced_parameters.dimensions = sliced_dimensions;
    sliced_parameters.bounds = sliced_bounds;
    sliced_parameters.smoothing_params = parameters.smoothing_params;
    sliced_parameters.flags = parameters.flags;
    sliced_parameters.instrument_type = parameters.instrument_type;

    for (size_t j = j_min; j <= j_max; ++j) {
        size_t offset = j * parameters.dimensions.n + i_min;
        stream.write(reinterpret_cast<const char *>(&source[offset]),
                     sizeof(double) * sliced_parameters.dimensions.n);
    }
    return Grid::Files::Dat::write_parameters(stream, sliced_parameters);
}

bool Grid::Serialize::read_parameters(std::istream &stream,
                                      Grid::Parameters *parameters) {
    Grid::Parameters read_params = {};
    Serialization::read_uint32(stream, &read_params.dimensions.n);
    Serialization::read_uint32(stream, &read_params.dimensions.m);
    Serialization::read_double(stream, &read_params.bounds.min_rt);
    Serialization::read_double(stream, &read_params.bounds.max_rt);
    Serialization::read_double(stream, &read_params.bounds.min_mz);
    Serialization::read_double(stream, &read_params.bounds.max_mz);
    Serialization::read_double(stream, &read_params.smoothing_params.mz);
    Serialization::read_double(stream, &read_params.smoothing_params.sigma_mz);
    Serialization::read_double(stream, &read_params.smoothing_params.sigma_rt);
    Serialization::read_uint8(stream, &read_params.instrument_type);
    Serialization::read_uint8(stream, &read_params.flags);
    *parameters = read_params;
    return stream.good();
}

bool Grid::Serialize::write_parameters(std::ostream &stream,
                                       const Grid::Parameters &parameters) {
    Serialization::write_uint32(stream, parameters.dimensions.n);
    Serialization::write_uint32(stream, parameters.dimensions.m);
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

bool Grid::Files::Dat::read_parameters(std::istream &stream,
                                       Grid::Parameters *parameters) {
    if (parameters == nullptr) {
        return false;
    }
    // Read the last byte of the stream to determine the length of the footer.
    stream.seekg(-1, std::ios::end);
    stream.seekg(-1 * stream.peek(), std::ios::end);
    // Read the parameters into the Grid::Parameters object. Note that we are
    // not making use of the Grid::Files::Dat::Parameters.spec_version yet, for
    // now we always read the data in the same way.
    Grid::Serialize::read_parameters(stream, parameters);
    return stream.good();
}

bool Grid::Files::Dat::write_parameters(std::ostream &stream,
                                        const Grid::Parameters &parameters) {
    Grid::Files::Dat::Parameters file_parameters = {
        0, Grid::Files::Dat::footer_size};
    Grid::Serialize::write_parameters(stream, parameters);
    Serialization::write_uint8(stream, file_parameters.spec_version);
    Serialization::write_uint8(stream, file_parameters.footer_size);
    return stream.good();
}

bool Grid::Files::Rawdump::write(std::ostream &stream,
                                 const std::vector<Grid::Point> &points) {
    uint64_t n_points = points.size();
    stream.write(reinterpret_cast<const char *>(&n_points), sizeof(uint64_t));
    stream.write(reinterpret_cast<const char *>(&points[0]),
                 sizeof(Grid::Point) * n_points);
    return stream.good();
}

bool Grid::Files::Rawdump::read(std::istream &stream,
                                std::vector<Grid::Point> &points) {
    uint64_t n_points = 0;
    stream.read(reinterpret_cast<char *>(&n_points), sizeof(uint64_t));
    if (stream.bad()) {
        return false;
    }
    points.resize(n_points);
    stream.read(reinterpret_cast<char *>(&points[0]),
                sizeof(Grid::Point) * n_points);
    return stream.good();
}
