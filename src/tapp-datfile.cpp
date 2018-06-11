#include "tapp-datfile.hpp"
#include <iostream>
#include <string>

// TODO: Should we have a function to verify if the parameters are correct?
bool DatFile::load(std::istream &stream, std::vector<double> *destination,
                   Grid::Parameters *parameters) {
    if (DatFile::load_parameters(stream, parameters)) {
        auto n_points = parameters->dimensions.n * parameters->dimensions.m;
        destination->resize(n_points);
        stream.seekg(0, std::ios::beg);
        stream.read(reinterpret_cast<char *>(&(*destination)[0]),
                    sizeof(double) * n_points);
    }
    return stream.good();
}

bool DatFile::write(std::ostream &stream, std::vector<double> &source,
                    Grid::Parameters &parameters) {
    stream.write(reinterpret_cast<char *>(&source[0]),
                 sizeof(double) * source.size());
    DatFile::write_parameters(stream, parameters);
    return stream.good();
}

bool DatFile::load_parameters(std::istream &stream,
                              Grid::Parameters *parameters) {
    // Read the last byte of the stream to determine the length of the footer.
    stream.seekg(-1, std::ios::end);
    stream.seekg(-1 * stream.peek(), std::ios::end);
    // Read the parameters into the Grid::Parameters object. Note that we are
    // not making use of the DatFile::Parameters.spec_version yet, for now we
    // always read the data in the same way.
    stream.read(reinterpret_cast<char *>(parameters), sizeof(Grid::Parameters));
    return stream.good();
}

bool DatFile::write_parameters(std::ostream &stream,
                               Grid::Parameters &parameters) {
    auto footer_size = static_cast<char>(sizeof(Grid::Parameters) +
                                         sizeof(DatFile::Parameters));
    DatFile::Parameters file_parameters = {1, footer_size};
    stream.write(reinterpret_cast<char *>(&parameters), sizeof(parameters));
    stream.write(reinterpret_cast<char *>(&file_parameters),
                 sizeof(file_parameters));
    return stream.good();
}

bool DatFile::load_uint32(std::istream &stream, uint32_t *i) {
    stream.read(reinterpret_cast<char *>(i), 4 * sizeof(char));
    return stream.good();
}

bool DatFile::save_uint32(std::ostream &stream, uint32_t i) {
    stream.write(reinterpret_cast<char *>(&i), 4 * sizeof(char));
    return stream.good();
}

bool DatFile::load_double(std::istream &stream, double *d) {
    stream.read(reinterpret_cast<char *>(d), sizeof(double));
    return stream.good();
}

bool DatFile::save_double(std::ostream &stream, double d) {
    stream.write(reinterpret_cast<char *>(&d), sizeof(double));
    return stream.good();
}
