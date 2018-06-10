#include "tapp-datfile.hpp"
#include <iostream>
#include <string>

bool DatFile::load(std::istream &stream, std::vector<double> &destination,
                   Grid::Parameters &parameters) {
    if (stream.good()) {
        // Reading all data from the stream into m_data.
        // m_data.resize(dimensions.n * dimensions.m);
        // m_dimensions = dimensions;
        // m_bounds = bounds;
        // m_instrument_type = instrument_type;
        // m_smoothing_params = smoothing_params;
        // if (stream.read(reinterpret_cast<char *>(&m_data[0]),
        // m_data.size() * sizeof(double))) {
        // return true;
        //}
    }
    return false;
}

bool DatFile::write(std::ostream &stream, const std::vector<double> &source,
                    const Grid::Parameters &parameters) {
    // if (stream.good()) {
    // stream.write(reinterpret_cast<char *>(&source[0]),
    // sizeof(double) * m_data.size());
    // return stream.good();
    //}
    return false;
}

struct FileFooter {
    std::string file_type = "TAPP";
    char footer_length;
};

bool DatFile::load_parameters(std::istream &stream,
                              Grid::Parameters *parameters) {
    if (stream.good()) {
        std::cout << "GOOD STREAM" << std::endl;
        std::cout << stream.tellg() << std::endl;
        stream.seekg(0, std::ios::beg);
        std::cout << stream.tellg() << std::endl;
        // auto stream_beg = stream.tellg();

        // Get to the end of the file to read the file footer.
        // stream.seekg(0, stream.end);
        // Read the length of the footer.
        char footer_length = 0;
        std::cout << "GOOD STREAM? " << stream.good() << std::endl;
        std::cout << "BEFORE: " << int(footer_length) << std::endl;
        footer_length = stream.peek();
        // stream.read(&footer_length, 1);
        std::cout << "AFTER: " << int(footer_length) << std::endl;
        std::cout << "GOOD STREAM? " << stream.good() << std::endl;
        std::cout << stream.tellg() << std::endl;
        // stream.seekg(0, stream.beg);
        std::cout << stream.tellg() << std::endl;
        return parameters;
    }
    return std::nullopt;
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
