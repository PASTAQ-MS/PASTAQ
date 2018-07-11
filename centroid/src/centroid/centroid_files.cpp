#include "centroid/centroid_files.hpp"

#include <cstring>

bool write_uint32(std::ostream &stream, uint32_t value) {
    // The data stream is encoded in little endian.
    uint8_t bytes[4];
    bytes[3] = value >> 24;
    bytes[2] = value >> 16;
    bytes[1] = value >> 8;
    bytes[0] = value >> 0;
    stream.write(reinterpret_cast<const char *>(&bytes), sizeof(uint32_t));
    return stream.good();
}

bool read_uint32(std::istream &stream, uint32_t *value) {
    // Read bytes from the stream.
    uint32_t read_value;
    stream.read(reinterpret_cast<char *>(&read_value), sizeof(uint32_t));

    // Interpret the bytes as little endian stream.
    uint8_t *bytes = reinterpret_cast<uint8_t *>(value);
    bytes[3] = read_value >> 24;
    bytes[2] = read_value >> 16;
    bytes[1] = read_value >> 8;
    bytes[0] = read_value >> 0;
    return stream.good();
}

bool write_uint64(std::ostream &stream, uint64_t value) {
    // The data stream is encoded in little endian.
    uint8_t bytes[8];
    bytes[7] = value >> 56;
    bytes[6] = value >> 48;
    bytes[5] = value >> 40;
    bytes[4] = value >> 32;
    bytes[3] = value >> 24;
    bytes[2] = value >> 16;
    bytes[1] = value >> 8;
    bytes[0] = value >> 0;
    stream.write(reinterpret_cast<const char *>(&bytes), sizeof(uint64_t));
    return stream.good();
}

bool read_uint64(std::istream &stream, uint64_t *value) {
    // Read bytes from the stream.
    uint64_t read_value;
    stream.read(reinterpret_cast<char *>(&read_value), sizeof(uint64_t));

    // Interpret the bytes as little endian stream.
    uint8_t *bytes = reinterpret_cast<uint8_t *>(value);
    bytes[7] = read_value >> 56;
    bytes[6] = read_value >> 48;
    bytes[5] = read_value >> 40;
    bytes[4] = read_value >> 32;
    bytes[3] = read_value >> 24;
    bytes[2] = read_value >> 16;
    bytes[1] = read_value >> 8;
    bytes[0] = read_value >> 0;
    return stream.good();
}

bool Centroid::Files::Bpks::write_double(std::ostream &stream, double value) {
    uint64_t raw_value;
    std::memcpy(&raw_value, &value, sizeof(value));
    return write_uint64(stream, raw_value);
}

bool Centroid::Files::Bpks::read_double(std::istream &stream, double *value) {
    stream.read(reinterpret_cast<char *>(value), 8);
    return stream.good();
    uint64_t raw_value;
    if (!read_uint64(stream, &raw_value)) {
        return false;
    }
    std::memcpy(value, &raw_value, sizeof(raw_value));
    return true;
}

bool Centroid::Files::Bpks::write_header(std::ostream &stream,
                                         const Header &header) {
    // auto footer_size = static_cast<char>(sizeof(Grid::Parameters) +
    // sizeof(Grid::Files::Dat::Parameters));
    // Grid::Files::Dat::Parameters file_parameters = {1, footer_size};
    // stream.write(reinterpret_cast<const char *>(&parameters),
    // sizeof(parameters));
    // stream.write(reinterpret_cast<const char *>(&file_parameters),
    // sizeof(file_parameters));
    return stream.good();
}
