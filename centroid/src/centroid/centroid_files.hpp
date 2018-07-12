#ifndef CENTROID_CENTROIDFILES_HPP
#define CENTROID_CENTROIDFILES_HPP
#include <iostream>

#include "centroid/centroid.hpp"

namespace Centroid::Files::Bpks {

// The header contains the necessary information for proper interpretation of
// the serialized data.
struct Header {
    uint8_t spec_version;
    uint64_t num_peaks;
    Grid::Parameters grid_params;
};

// Read/Write the header from the given binary stream.
bool read_header(std::istream &stream, Header *header);
bool write_header(std::ostream &stream, const Header &header);

// TODO(alex): Move to Centroid::Serialization?
// Read/Write a single peak from the given binary stream.
bool read_peak(std::istream &stream, Centroid::Peak &peak);
bool write_peak(std::ostream &stream, const Centroid::Peak &peak);

// TODO(alex): Move to Centroid::Serialization?
// Read/Write all peaks from the given binary stream.
bool read_peaks(std::istream &stream, Header &header,
                std::vector<Centroid::Peak> &peak);
bool write_peaks(std::ostream &stream, Header &header,
                 const std::vector<Centroid::Peak> &peak);

}  // namespace Centroid::Files::Bpks

namespace Centroid::Files::Csv {
// TODO(alex): Dump peaks to csv.
// TODO(alex): Read peaks from csv.
}  // namespace Centroid::Files::Csv

#endif /* CENTROID_CENTROIDFILES_HPP */
