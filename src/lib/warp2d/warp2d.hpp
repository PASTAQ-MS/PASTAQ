#ifndef WARP2D_WARP2D_HPP
#define WARP2D_WARP2D_HPP

#include <cstdint>
#include <vector>

#include "centroid/centroid.hpp"

namespace Warp2D {
// The parameters used in Warp2D.
//
// - slack: Determines how much each node can be warped in either direction.
// - window_size: The number of points for each segment.
// - num_points: The number of points in which the retention time range will be
//   divided.
// - peaks_per_window: The number of peaks that will be used in each of the
//   segments.
// - rt_expand_factor: The algorithm will not warp at the range limits so we
//   need to expand the range by some factor.
// TODO: Revise if the type of these parameters is correct, i.e. Why would the
// slack be negative?
struct Parameters {
    int64_t slack;
    int64_t window_size;
    int64_t num_points;
    int64_t peaks_per_window;
    double rt_expand_factor;
};

// The main element of the FU matrix used by the Correlation Optimised Warping
// (COW) algorithm. In the original algorithm it is described the use of two
// matrices F and U to store the cumulative correlation/similarity and the
// optimal precursor for a given element, but since the algorithm doesn't
// require the majority of the elements of that matrix it makes more sense to
// structure the data as a list FU nodes for each segment.
struct Node {
    double f;   // Cumulative similarity value.
    int64_t u;  // Optimal predecessor index.
};

// This structure describes the required data to perform a warping between two
// nodes.
struct PotentialWarping {
    int64_t i;          // Index of the node on the current level for x_start.
    int64_t j;          // Index of the node on the next level for x_end.
    int64_t src_start;  // The initial point to warp.
    int64_t src_end;    // The end point to warp.
    double warped_similarity;  // The similarity obtained after warping.
};

// Each level contains the potential warpings for the segment, the FU nodes for
// optimal warping detection and the start and end points of the level. For
// example, given the following parameters:
//
//     slack = 1
//     window_size = 5
//     num_points = 25
//
// We will have 5 segments with the nodes distributed in the following way (Note
// that we will always have 1 extra segment at the end:
//
//              0         1         2         3         4
//     Level 0 |x| | | | | | | | | | | | | | | | | | | | | | | | |
//     Level 1 | | | | |x|x|x| | | | | | | | | | | | | | | | | | |
//     Level 2 | | | | | | | | |x|x|x|x|x| | | | | | | | | | | | |
//     Level 3 | | | | | | | | | | | | | |x|x|x|x|x| | | | | | | |
//     Level 4 | | | | | | | | | | | | | | | | | | | |x|x|x| | | |
//     Level 5 | | | | | | | | | | | | | | | | | | | | | | | | |x|
//
// - Level 0:
//     - start: 0
//     - end: 0
//     - potential_warpings:
//         [
//             i: 0 j: 0 src_start: 0 src_end: 4
//             i: 0 j: 1 src_start: 0 src_end: 5
//             i: 0 j: 2 src_start: 0 src_end: 6
//         ]
// - level: 1
//     - start: 4
//     - end: 6
//     - potential_warpings:
//         [
//              i: 0 j: 0 src_start: 4 src_end: 8
//              i: 0 j: 1 src_start: 4 src_end: 9
//              i: 0 j: 2 src_start: 4 src_end: 10
//              i: 1 j: 1 src_start: 5 src_end: 9
//              i: 1 j: 2 src_start: 5 src_end: 10
//              i: 1 j: 3 src_start: 5 src_end: 11
//              i: 2 j: 2 src_start: 6 src_end: 10
//              i: 2 j: 3 src_start: 6 src_end: 11
//              i: 2 j: 4 src_start: 6 src_end: 12
//         ]
// - level: 2
//     - start: 8
//     - end: 12
//     - potential_warpings:
//         [
//              i: 0 j: 0 src_start: 8 src_end: 13
//              i: 0 j: 1 src_start: 8 src_end: 14
//              i: 1 j: 0 src_start: 9 src_end: 13
//              i: 1 j: 1 src_start: 9 src_end: 14
//              i: 1 j: 2 src_start: 9 src_end: 15
//              i: 2 j: 1 src_start: 10 src_end: 14
//              i: 2 j: 2 src_start: 10 src_end: 15
//              i: 2 j: 3 src_start: 10 src_end: 16
//              i: 3 j: 2 src_start: 11 src_end: 15
//              i: 3 j: 3 src_start: 11 src_end: 16
//              i: 3 j: 4 src_start: 11 src_end: 17
//              i: 4 j: 3 src_start: 12 src_end: 16
//              i: 4 j: 4 src_start: 12 src_end: 17
//         ]
// - level: 3
//     - start: 13
//     - end: 17
//     - potential_warpings:
//         [
//              i: 0 j: 0 src_start: 13 src_end: 19
//              i: 1 j: 0 src_start: 14 src_end: 19
//              i: 1 j: 1 src_start: 14 src_end: 20
//              i: 2 j: 0 src_start: 15 src_end: 19
//              i: 2 j: 1 src_start: 15 src_end: 20
//              i: 2 j: 2 src_start: 15 src_end: 21
//              i: 3 j: 1 src_start: 16 src_end: 20
//              i: 3 j: 2 src_start: 16 src_end: 21
//              i: 4 j: 2 src_start: 17 src_end: 21
//         ]
// - level: 4
//     - start: 19
//     - end: 21
//     - potential_warpings:
//         [
//              i: 0 j: 0 src_start: 19 src_end: 25
//              i: 1 j: 0 src_start: 20 src_end: 25
//              i: 2 j: 0 src_start: 21 src_end: 25
//         ]
// - level: 5
//     - start: 25
//     - end: 25
//     - potential_warpings: []
//
struct Level {
    int64_t start;
    int64_t end;
    std::vector<Node> nodes;
    std::vector<PotentialWarping> potential_warpings;
};

struct TimeMap {
    uint64_t num_segments;
    double rt_min;
    double rt_max;
    std::vector<double> rt_start;
    std::vector<double> rt_end;
    std::vector<double> sample_rt_start;
    std::vector<double> sample_rt_end;
};

// Warp the peaks by linearly interpolating their retention time to the given
// reference time. Note that we are just performing linear displacement of the
// center of the peaks, we do not deform the peak shape by adjusting the sigmas.
std::vector<Centroid::Peak> interpolate_peaks(
    const std::vector<Centroid::Peak>& source_peaks, double source_rt_start,
    double source_rt_end, double ref_rt_start, double ref_rt_end);

// Returns a copy of the peaks from source_peaks that are in the given region
// between time_start and time_end.
std::vector<Centroid::Peak> peaks_in_rt_range(
    const std::vector<Centroid::Peak>& source_peaks, double time_start,
    double time_end);

// Filter the peaks based on peak height. Note that this function modifies the
// given `peaks` argument by sorting the vector in place.
std::vector<Centroid::Peak> filter_peaks(std::vector<Centroid::Peak>& peaks,
                                         size_t n_peaks_max);

// Initialize the vector of Levels, including the potential warpings and FU
// nodes.
std::vector<Level> initialize_levels(int64_t num_sectors, int64_t window_size,
                                     int64_t slack, int64_t num_points);

// Calculate all warped similarities from each PotentialWarping in
// level.warped_similarities.
void compute_warped_similarities(
    Warp2D::Level& level, double rt_start, double rt_end, double rt_min,
    double delta_rt, const std::vector<Centroid::Peak>& ref_peaks,
    const std::vector<Centroid::Peak>& source_peaks);

// Calculates the optimal set of warping points using the computed warped
// similarities in levels. It does so in two steps: First it walks back the list
// of warped similarities and updates the FU nodes, and then it walks forward
// the FU nodes to find the optimal warping path.
std::vector<int64_t> find_optimal_warping(std::vector<Level>& levels);

// Perform the Warp2D algorithm to find the optimal TimeMap for peak warping.
TimeMap calculate_time_map(const std::vector<Centroid::Peak>& ref_peaks,
                           const std::vector<Centroid::Peak>& source_peaks,
                           const Parameters& parameters, uint64_t max_threads);

// Use the given TimeMap to interpolate the source_peaks for retention time
// alignment.
std::vector<Centroid::Peak> warp_peaks(
    const std::vector<Centroid::Peak>& source_peaks, const TimeMap& time_map);

// Use a TimeMap to interpolate a given retention time.
double warp(const TimeMap& time_map, double rt);

}  // namespace Warp2D

#endif /* WARP2D_WARP2D_HPP */
