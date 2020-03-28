#ifndef FEATUREDETECTION__FEATUREDETECTION_HPP
#define FEATUREDETECTION__FEATUREDETECTION_HPP
#include <optional>

#include "link/link.hpp"
#include "utils/search.hpp"

namespace FeatureDetection {
// TODO: Expand documentation.
struct Feature {
    uint64_t id;
    double score;
    double average_rt;
    double average_rt_delta;
    double average_rt_sigma;
    double average_mz;  // NOTE: This is a weighted average, maybe we should
                        // change the name to be more explicit.
    double average_mz_sigma;
    double total_height;
    double total_volume;
    double max_height;
    double max_volume;
    double monoisotopic_mz;
    double monoisotopic_rt;
    double monoisotopic_height;
    double monoisotopic_volume;
    int8_t charge_state;  // FIXME: Why int8 instead of uint8?
    std::vector<uint64_t> peak_ids;
};

// TODO: Expand documentation.
struct TheoreticalIsotopes {
    std::vector<double> mzs;
    std::vector<double> percs;
};

struct RootNode {
    std::vector<uint64_t> nodes_next;
    std::vector<uint64_t> nodes_prev;
    bool visited;
    uint64_t id;
};

typedef std::vector<RootNode> CandidateGraph;

std::vector<Feature> detect_features(const std::vector<Centroid::Peak> &peaks,
                                     const std::vector<uint8_t> &charge_states);

}  // namespace FeatureDetection

#endif /* FEATUREDETECTION__FEATUREDETECTION_HPP */
