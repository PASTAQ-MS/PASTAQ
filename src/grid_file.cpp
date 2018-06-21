#include <iostream>
#include <string>

#include "grid_file.hpp"

bool Grid::File::load(std::istream &stream, std::vector<double> *destination,
                      Grid::Parameters *parameters) {
    if (Grid::File::load_parameters(stream, parameters)) {
        auto n_points = parameters->dimensions.n * parameters->dimensions.m;
        destination->resize(n_points);
        stream.seekg(0, std::ios::beg);
        stream.read(reinterpret_cast<char *>(&(*destination)[0]),
                    sizeof(double) * n_points);
    }
    return stream.good();
}

bool Grid::File::write(std::ostream &stream, const std::vector<double> &source,
                       const Grid::Parameters &parameters) {
    stream.write(reinterpret_cast<const char *>(&source[0]),
                 sizeof(double) * source.size());
    Grid::File::write_parameters(stream, parameters);
    return stream.good();
}

bool Grid::File::load_parameters(std::istream &stream,
                                 Grid::Parameters *parameters) {
    // Read the last byte of the stream to determine the length of the footer.
    stream.seekg(-1, std::ios::end);
    stream.seekg(-1 * stream.peek(), std::ios::end);
    // Read the parameters into the Grid::Parameters object. Note that we are
    // not making use of the Grid::File::Parameters.spec_version yet, for now we
    // always read the data in the same way.
    stream.read(reinterpret_cast<char *>(parameters), sizeof(Grid::Parameters));
    return stream.good();
}

bool Grid::File::write_parameters(std::ostream &stream,
                                  const Grid::Parameters &parameters) {
    auto footer_size = static_cast<char>(sizeof(Grid::Parameters) +
                                         sizeof(Grid::File::Parameters));
    Grid::File::Parameters file_parameters = {1, footer_size};
    stream.write(reinterpret_cast<const char *>(&parameters), sizeof(parameters));
    stream.write(reinterpret_cast<const char *>(&file_parameters),
                 sizeof(file_parameters));
    return stream.good();
}
