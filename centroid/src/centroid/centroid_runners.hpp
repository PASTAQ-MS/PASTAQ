#ifndef CENTROID_CENTROIDRUNNERS_HPP
#define CENTROID_CENTROIDRUNNERS_HPP

#include <vector>

#include "centroid.hpp"

namespace Centroid::Runners::Serial {

// Find the peaks in serial.
std::vector<Centroid::Peak> run(const RawData::RawData &raw_data,
                                const Grid::Mesh &mesh, size_t max_peaks);

}  // namespace Centroid::Runners::Serial

namespace Centroid::Runners::Parallel {

// Find the peaks in parallel.
std::vector<Centroid::Peak> run(const RawData::RawData &raw_data,
                                const Grid::Mesh &mesh, size_t max_peaks);

}  // namespace Centroid::Runners::Parallel

#endif /* CENTROID_CENTROIDRUNNERS_HPP */
