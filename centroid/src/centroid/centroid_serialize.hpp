#ifndef CENTROID_CENTROIDSERIALIZE_HPP
#define CENTROID_CENTROIDSERIALIZE_HPP

#include <iostream>

#include "centroid/centroid.hpp"

// This namespace groups the functions used to serialize Centroid data
// structures into a binary stream.
namespace Centroid::Serialize {

// Centroid::Peak
bool read_peak(std::istream &stream, Centroid::Peak *peak);
bool write_peak(std::ostream &stream, const Centroid::Peak &peak);

}  // namespace Centroid::Serialize

#endif /* CENTROID_CENTROIDSERIALIZE_HPP */
