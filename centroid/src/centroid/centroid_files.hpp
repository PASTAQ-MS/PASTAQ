#ifndef CENTROID_CENTROIDFILES_HPP
#define CENTROID_CENTROIDFILES_HPP
#include <iostream>

#include "centroid/centroid.hpp"

namespace Centroid::Files::Bpks {

struct Header {
    uint64_t header_length;
    uint64_t num_peaks;
    unsigned char spec_version;
    Grid::Parameters grid_parameters;
};

// Read/Write the header from the given binary stream.
bool read_header(std::istream &stream, Header &header);
bool write_header(std::ostream &stream, const Header &header);

// Read/Write a single peak from the given binary stream.
bool read_peak(std::istream &stream, Centroid::Peak &peak);
bool write_peak(std::ostream &stream, const Centroid::Peak &peak);

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
