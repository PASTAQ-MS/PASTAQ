#ifndef CENTROID_CENTROIDSERIALIZE_HPP
#define CENTROID_CENTROIDSERIALIZE_HPP

#include <iostream>

#include "centroid/centroid.hpp"

// This namespace groups the functions used to serialize Centroid data
// structures into a binary stream.
namespace Centroid::Serialize {

// Read/Write a single peak to/from the given binary stream.
bool read_peak(std::istream &stream, Centroid::Peak *peak);
bool write_peak(std::ostream &stream, const Centroid::Peak &peak);

// Read/Write all peaks to/from the given binary stream.
bool read_peaks(std::istream &stream, std::vector<Centroid::Peak> *peaks);
bool write_peaks(std::ostream &stream,
                 const std::vector<Centroid::Peak> &peaks);

}  // namespace Centroid::Serialize

#endif /* CENTROID_CENTROIDSERIALIZE_HPP */
