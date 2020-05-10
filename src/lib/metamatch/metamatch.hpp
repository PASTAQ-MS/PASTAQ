#ifndef METAMATCH_METAMATCH_HPP
#define METAMATCH_METAMATCH_HPP

#include <cstdint>
#include <vector>

#include "centroid/centroid.hpp"
#include "feature_detection/feature_detection.hpp"

namespace MetaMatch {

// In order for a cluster to be considered valid a number of peaks on any given
// class must appear on the cluster. For example if we have 10 files in class
// `A` and 20 in class `B` and the required fraction is `0.5`, a cluster will be
// considered valid if it has peaks coming from 5 class `A` files or 10 class
// `B` files.
struct ClassMap {
    uint32_t id;
    uint32_t n_files;
    uint32_t required_hits;
};

// A MetaMatch::Peak is an extension of Centroid::Peak that allow us to store
// the necessary information for the clustering algorithm.
// FIXME: This might be overkill, I must find a better way of doing this.
// Essentially we are requiring the copy of ALL the peaks in ALL the files for
// annotation. It might not be a big deal with the file sizes we are currently
// using, but why this waste?
struct Peak : Centroid::Peak {
    uint32_t file_id;
    uint32_t class_id;
    // FIXME: We are sacrificing half of our address space just so that we can
    // put a -1 if the cluster is not linked yet, perhaps it's better to keep
    // track of this as a separate vector so that we can use uint64_t
    int64_t cluster_id;
    double cluster_mz;
    double cluster_rt;
};

// The aggregate information for the peaks of a given cluster.
struct FeatureId {
    uint64_t file_id;
    uint64_t feature_id;
};
// TODO: Make this general in order for it to work with features OR peaks.
struct Cluster {
    // FIXME: We are sacrificing half of our address space just so that we can
    // put a -1 if the cluster is not linked yet, perhaps it's better to keep
    // track of this as a separate vector so that we can use uint64_t
    int64_t id;
    double mz;
    double rt;
    // TODO(alex): Other stats here...
    std::vector<double> file_heights;
    std::vector<double> file_volumes;
    double avg_height;

    // The feature ids on each file associated with this cluster.
    // FIXME: Naming peak_ids/feature_ids!
    std::vector<FeatureId> peak_ids;
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

// Performs a greedy feature clustering algorithm by trying to match features
// from multiple files based on their charge state and monoisotopic peak. The
// features are searched in descending intensity order.
std::vector<FeatureCluster> find_feature_clusters(
    std::vector<uint64_t>& group_ids,
    std::vector<std::vector<FeatureDetection::Feature>>& features,
    double keep_perc, double intensity_threshold, double n_sig_mz,
    double n_sig_rt);

// Performs a centroid based clustering algorithm. This algorithm modifies the
// given peaks array in place, changing the order and the clustering
// information.
void find_clusters(std::vector<MetaMatch::Peak>& peaks,
                   const std::vector<ClassMap>& class_maps,
                   double n_sig_mz, double n_sig_rt);

// Returns the orphan peaks and remove them from the given peaks array.
std::vector<MetaMatch::Peak> extract_orphans(
    std::vector<MetaMatch::Peak>& peaks);

// Gather the peaks belonging to a cluster into a MetaMatch::Cluster. Since the
// peaks need to be sorted prior to perform the aggregation, the original peaks
// array order might change.
std::vector<MetaMatch::Cluster> reduce_cluster(
    std::vector<MetaMatch::Peak>& peaks, size_t n_files);

}  // namespace MetaMatch

#endif /* METAMATCH_METAMATCH_HPP */
