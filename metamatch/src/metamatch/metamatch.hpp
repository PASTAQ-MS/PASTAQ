#ifndef METAMATCH_METAMATCH_HPP
#define METAMATCH_METAMATCH_HPP

#include <cstdint>
#include <vector>

#include "centroid/centroid_files.hpp"

namespace MetaMatch {

// TODO(alex): Add docs...
struct Peak : Centroid::Peak {
    size_t file_id;
    size_t class_id;
    int cluster_id;
    double cluster_mz;
    double cluster_rt;
};

struct Cluster {
    int id;
    double mz;
    double rt;
    // TODO(alex): Other stats here...
    std::vector<double> file_heights;
};

// TODO(alex): Add docs...
void find_candidates(std::vector<MetaMatch::Peak>& peaks, double radius_mz,
                     double radius_rt, size_t n_files, size_t n_classes,
                     double fraction);

// TODO(alex): Add docs...
std::vector<MetaMatch::Peak> extract_orphans(
    std::vector<MetaMatch::Peak>& peaks);
std::vector<MetaMatch::Cluster> reduce_cluster(
    std::vector<MetaMatch::Peak>& peaks, size_t n_files);
}  // namespace MetaMatch

#endif /* METAMATCH_METAMATCH_HPP */
