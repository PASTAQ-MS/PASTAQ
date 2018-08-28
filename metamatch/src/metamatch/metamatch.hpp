#ifndef METAMATCH_METAMATCH_HPP
#define METAMATCH_METAMATCH_HPP

#include <cstdint>
#include <vector>

#include "centroid/centroid_files.hpp"

namespace MetaMatch {

struct Peak : Centroid::Peak {
    size_t file_id;
    size_t class_id;
    int cluster_id;
    double cluster_mz;
    double cluster_rt;
};

void find_candidates(std::vector<MetaMatch::Peak>& peaks, double radius_mz,
                     double radius_rt, size_t n_files, size_t n_classes,
                     double fraction);

std::vector<MetaMatch::Peak> extract_orphans(std::vector<MetaMatch::Peak>& peaks);
}  // namespace MetaMatch

#endif /* METAMATCH_METAMATCH_HPP */
