#ifndef WARP2D_WARP2DRUNNERS_HPP
#define WARP2D_WARP2DRUNNERS_HPP

#include <vector>

#include "warp2d/warp2d.hpp"

namespace Warp2D::Runners::Serial {

// Find the peaks on the given Grid in serial.
std::vector<Centroid::Peak> run(
    const std::vector<Centroid::Peak>& reference_peaks,
    const std::vector<Centroid::Peak>& source_peaks,
    const Warp2D::Parameters& parameters);

}  // namespace Warp2D::Runners::Serial

namespace Warp2D::Runners::Parallel {

// Find the peaks on the given Grid in parallel.
std::vector<Centroid::Peak> run(
    const std::vector<Centroid::Peak>& reference_peaks,
    const std::vector<Centroid::Peak>& source_peaks,
    const Warp2D::Parameters& parameters, uint64_t max_threads);

}  // namespace Warp2D::Runners::Parallel

#endif /* WARP2D_WARP2DRUNNERS_HPP */
