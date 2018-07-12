#include "centroid/centroid_files.hpp"
#include "centroid/centroid_serialize.hpp"
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

bool Centroid::Files::Bpks::read_peaks(std::istream &stream,
                                       Grid::Parameters *grid_parameters,
                                       std::vector<Centroid::Peak> *peaks) {
    Header header = {};
    if (!read_header(stream, &header)) {
        return false;
    }
    *grid_parameters = header.grid_params;
    peaks->resize(header.num_peaks);
    for (auto &peak : *peaks) {
        if (!Centroid::Serialize::read_peak(stream, &peak)) {
            return false;
        }
    }
    return stream.good();
}

bool Centroid::Files::Bpks::write_peaks(
    std::ostream &stream, const Grid::Parameters &grid_parameters,
    const std::vector<Centroid::Peak> &peaks) {
    if (!write_header(stream, {0, peaks.size(), grid_parameters})) {
        return false;
    }
    for (const auto &peak : peaks) {
        if (!Centroid::Serialize::write_peak(stream, peak)) {
            return false;
        }
    }
    return stream.good();
}
