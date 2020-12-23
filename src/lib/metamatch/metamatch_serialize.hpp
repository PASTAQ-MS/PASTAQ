#ifndef METAMATCH_METAMATCHSERIALIZE_HPP
#define METAMATCH_METAMATCHSERIALIZE_HPP

#include <iostream>

#include "metamatch/metamatch.hpp"

// This namespace groups the functions used to serialize MetaMatch data
// structures into a binary stream.
namespace MetaMatch::Serialize {

// MetaMatch::PeakCluster
bool read_peak_clusters(std::istream &stream,
                           std::vector<PeakCluster> *clusters);
bool write_peak_clusters(std::ostream &stream,
                        const std::vector<PeakCluster> &clusters);

// MetaMatch::FeatureCluster
bool read_feature_clusters(std::istream &stream,
                           std::vector<FeatureCluster> *clusters);
bool write_feature_clusters(std::ostream &stream,
                            const std::vector<FeatureCluster> &clusters);

}  // namespace MetaMatch::Serialize

#endif /* METAMATCH_METAMATCHSERIALIZE_HPP */
