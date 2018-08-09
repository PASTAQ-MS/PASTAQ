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
    int peaks_per_window;
};

// TODO(alex): add more docs.
struct Node {
    double f;  // Cummulative similarity value.
    int u;     // Optimal predecessor index.
};

// TODO(alex): add more docs.
struct PotentialWarping {
    int i;  // Index of the node on the current level for x_start.
    int j;  // Index of the node on the next level for x_end.
    int x_start;
    int x_end;
    double warped_similarity;
};

// TODO(alex): add more docs.
struct Level {
    int start;
    int end;
    std::vector<Node> nodes;
    std::vector<PotentialWarping> potential_warpings;
};

// Calculate the overlaping area between two peaks.
double peak_overlap(const Centroid::Peak& peak_a, const Centroid::Peak& peak_b);

// Calculate the cummulative similarity between two sets of peaks.
double similarity_2D(const std::vector<Centroid::Peak>& set_a,
                     const std::vector<Centroid::Peak>& set_b);

// FIXME: COMMENT NOT CORRECT... Warp the sample_peaks to target_peaks in the
// retention time dimension. The warping is performed by using a variant of the
// Correlation Optimised Warping (COW) that uses the overlaping volume of the
// peaks as the similarity/benefit function. Returns the peaks after successful
// warping.
std::vector<Centroid::Peak> warp_peaks(
    const std::vector<Centroid::Peak>& source_peaks, double source_rt_start,
    double source_rt_end, double ref_rt_start, double ref_rt_end);

// TODO(alex): add more docs.
std::vector<Centroid::Peak> peaks_in_rt_range(
    const std::vector<Centroid::Peak>& source_peaks, double time_start,
    double time_end);

// TODO(alex): add more docs.
std::vector<Centroid::Peak> filter_peaks(std::vector<Centroid::Peak>& peaks,
                                         size_t n_peaks_max);

// TODO(alex): add more docs.
std::vector<Level> initialize_levels(int N, int m, int t, int nP);

// TODO(alex): add more docs.
void compute_warped_similarities(
    Warp2D::Level& level, double rt_start, double rt_end, double rt_min,
    double delta_rt, const std::vector<Centroid::Peak>& target_peaks,
    const std::vector<Centroid::Peak>& source_peaks);

// TODO(alex): add more docs.
std::vector<int> find_optimal_warping(std::vector<Level>& levels, int N);
}  // namespace Warp2D

#endif /* WARP2D_WARP2D_HPP */
