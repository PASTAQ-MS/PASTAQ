#ifndef METAMATCH_METAMATCH_HPP
#define METAMATCH_METAMATCH_HPP

#include <cstdint>
#include <vector>

#include "centroid/centroid_files.hpp"

namespace MetaMatch {
struct Index {
    size_t file_index;
    size_t peak_index;
};

struct Peak : Centroid::Peak {
    size_t file_id;
    size_t class_id;
    int cluster_id;
    double cluster_mz;
    double cluster_rt;
};

void find_candidates(std::vector<MetaMatch::Peak>& peak_files, double radius_mz,
                     double radius_rt);
}  // namespace MetaMatch

#endif /* METAMATCH_METAMATCH_HPP */
