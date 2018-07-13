#ifndef CENTROID_CENTROIDRUNNERS_HPP
#define CENTROID_CENTROIDRUNNERS_HPP

#include <vector>

#include "centroid.hpp"

namespace Centroid::Runners::Serial {

// Find the peaks on the given Grid in serial.
std::vector<Centroid::Peak> run(const Centroid::Parameters &parameters,
                                const std::vector<double> &data);

}  // namespace Centroid::Runners::Serial

namespace Centroid::Runners::Parallel {

// Find the peaks on the given Grid in parallel.
std::vector<Centroid::Peak> run(uint64_t max_threads,
                                const Centroid::Parameters &parameters,
                                const std::vector<double> &data);

}  // namespace Centroid::Runners::Parallel

#endif /* CENTROID_CENTROIDRUNNERS_HPP */
