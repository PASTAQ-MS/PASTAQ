#ifndef METAMATCH_METAMATCH_HPP
#define METAMATCH_METAMATCH_HPP

#include <cstdint>
#include <vector>

#include "centroid/centroid.hpp"
#include "feature_detection/feature_detection.hpp"

namespace MetaMatch {

// The aggregate information for the peaks of a given cluster.
struct PeakId {
    uint64_t file_id;
    uint64_t peak_id;
};
struct FeatureId {
    uint64_t file_id;
    uint64_t feature_id;
};

struct PeakCluster {
    uint64_t id;
    double mz;
    double rt;
    // Statistics for this cluster.
    double avg_height;
    double avg_volume;
    // Quantifications for this cluster.
    std::vector<double> heights;
    std::vector<double> volumes;
    // The peak ids on each file associated with this cluster.
    std::vector<PeakId> peak_ids;
};

struct FeatureCluster {
    uint64_t id;
    double mz;
    double rt;
    int8_t charge_state;
    // Statistics for this cluster.
    double avg_total_height;
    double avg_monoisotopic_height;
    double avg_max_height;
    double avg_total_volume;
    double avg_monoisotopic_volume;
    double avg_max_volume;
    // Quantifications for this cluster.
    std::vector<double> total_heights;
    std::vector<double> monoisotopic_heights;
    std::vector<double> max_heights;
    std::vector<double> total_volumes;
    std::vector<double> monoisotopic_volumes;
    std::vector<double> max_volumes;
    // The feature ids on each file associated with this cluster.
    std::vector<FeatureId> feature_ids;
};

// Performs a greedy feature clustering algorithm by trying to match peaks or
// features from multiple files based on their charge state and monoisotopic
// peak (for the later). The peaks/features are searched in descending intensity
// order.
std::vector<MetaMatch::PeakCluster> find_peak_clusters(
    std::vector<uint64_t>& group_ids,
    std::vector<std::vector<Centroid::Peak>>& peaks,
    double keep_perc, double intensity_threshold, double n_sig_mz,
    double n_sig_rt);
std::vector<FeatureCluster> find_feature_clusters(
    std::vector<uint64_t>& group_ids,
    std::vector<std::vector<FeatureDetection::Feature>>& features,
    double keep_perc, double intensity_threshold, double n_sig_mz,
    double n_sig_rt);

}  // namespace MetaMatch

#endif /* METAMATCH_METAMATCH_HPP */
