#ifndef WARP2D_WARP2D_HPP
#define WARP2D_WARP2D_HPP

#include <cstdint>
#include <vector>

#include "centroid/centroid.hpp"

namespace Warp2D {
// FIXME(alex): Reference function. To be removed.
int cow_2D(std::vector<Centroid::Peak>& target_peaks,
           std::vector<Centroid::Peak>& source_peaks, int sample_length,
           int segment_length, int slack);

// Warp the sample_peaks to target_peaks in the retention time dimension. The
// warping is performed by using a variant of the Correlation Optimised Warping
// (COW) that uses the overlaping volume of the peaks as the similarity/benefit
// function.
int warp_peaks(std::vector<Centroid::Peak>& target_peaks,
               std::vector<Centroid::Peak>& source_peaks, int sample_length,
               int segment_length, int slack);
}  // namespace Warp2D

#endif /* WARP2D_WARP2D_HPP */
