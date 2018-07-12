#include "centroid/centroid_files.hpp"
#include "grid/grid_serialize.hpp"
#include "utils/serialization.hpp"

bool Centroid::Files::Bpks::write_header(std::ostream &stream,
                                         const Header &header) {
    Serialization::write_uint8(stream, header.spec_version);
    Serialization::write_uint64(stream, header.num_peaks);
    Grid::Serialize::write_parameters(stream, header.grid_params);
    return stream.good();
}

bool Centroid::Files::Bpks::read_header(std::istream &stream, Header *header) {
    Serialization::read_uint8(stream, &header->spec_version);
    Serialization::read_uint64(stream, &header->num_peaks);
    Grid::Serialize::read_parameters(stream, &header->grid_params);
    return stream.good();
}
