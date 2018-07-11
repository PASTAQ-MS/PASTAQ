#include "centroid/centroid_files.hpp"
#include "utils/serialization.hpp"

bool Centroid::Files::Bpks::write_header(std::ostream &stream,
                                         const Header &header) {
    Serialization::write_uint64(stream, header.header_length);
    Serialization::write_uint64(stream, header.num_peaks);
    Serialization::write_uint8(stream, header.spec_version);
    // TODO(alex): Write Grid::Parameters
    return stream.good();
}

bool Centroid::Files::Bpks::read_header(std::istream &stream, Header *header) {
    Header read_header = {};
    Serialization::read_uint64(stream, &read_header.header_length);
    Serialization::read_uint64(stream, &read_header.num_peaks);
    Serialization::read_uint8(stream, &read_header.spec_version);
    // TODO(alex): Read Grid::Parameters
    *header = read_header;
    return stream.good();
}
