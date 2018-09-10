#ifndef METAMATCH_METAMATCH_HPP
#define METAMATCH_METAMATCH_HPP

#include <cstdint>
#include <vector>

#include "centroid/centroid_files.hpp"

namespace MetaMatch {

// In order for a cluster to be considered valid a number of peaks on any given
// class must appear on the cluster. For example if we have 10 files in class
// `A` and 20 in class `B` and the required fraction is `0.5`, a cluster will be
// considered valid if it has peaks coming from 5 class `A` files or 10 class
// `B` files.
struct ClassMap {
    size_t id;
    size_t n_files;
    size_t required_hits;
};

// TODO(alex): Currently we don't use adaptative mz radius.
struct Parameters {
    double radius_mz;
    double radius_rt;
    std::vector<ClassMap> class_maps;
};

// A MetaMatch::Peak is an extension of Centroid::Peak that allow us to store
// the necessary information for the clustering algorithm.
struct Peak : Centroid::Peak {
    size_t file_id;
    size_t class_id;
    int cluster_id;
    double cluster_mz;
    double cluster_rt;
};

// The aggregate information for the peaks of a given cluster.
struct Cluster {
    int id;
    double mz;
    double rt;
    // TODO(alex): Other stats here...
    std::vector<double> file_heights;
};

// Performs a centroid based clustering algorithm. This algorithm modifies the
// given peaks array in place, changing the order and the clustering
// information.
void find_clusters(std::vector<MetaMatch::Peak>& peaks,
                   const MetaMatch::Parameters& parameters);

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
