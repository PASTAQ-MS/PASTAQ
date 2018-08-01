#ifndef WARP2D_WARP2D_HPP
#define WARP2D_WARP2D_HPP

#include <cstdint>
#include <vector>

#include "centroid/centroid.hpp"

namespace Warp2D {
// The parameters used in Warp2D
// TODO(alex): add more docs.
struct Parameters {
    int slack;
    int window_size;
    int num_points;
};

// Calculate the overlaping area between two peaks.
double peak_overlap(const Centroid::Peak& peak_a, const Centroid::Peak& peak_b);

// Calculate the cummulative similarity between two sets of peaks.
double similarity_2D(const std::vector<Centroid::Peak>& set_a,
                     const std::vector<Centroid::Peak>& set_b);

// Warp the sample_peaks to target_peaks in the retention time dimension. The
// warping is performed by using a variant of the Correlation Optimised Warping
// (COW) that uses the overlaping volume of the peaks as the similarity/benefit
// function. Returns the peaks after successful warping.
std::vector<Centroid::Peak> warp_peaks(
    const std::vector<Centroid::Peak>& target_peaks,
    const std::vector<Centroid::Peak>& source_peaks,
    const Warp2D::Parameters& parameters);
}  // namespace Warp2D

#endif /* WARP2D_WARP2D_HPP */
