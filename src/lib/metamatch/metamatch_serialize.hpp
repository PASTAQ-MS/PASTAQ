#ifndef METAMATCH_METAMATCHSERIALIZE_HPP
#define METAMATCH_METAMATCHSERIALIZE_HPP

#include <iostream>

#include "metamatch/metamatch.hpp"

// This namespace groups the functions used to serialize MetaMatch data
// structures into a binary stream.
namespace MetaMatch::Serialize {

// MetaMatch::Cluster
bool read_cluster(std::istream &stream, Cluster *cluster);
bool write_cluster(std::ostream &stream, const Cluster &cluster);

// std::vector<MetaMatch::Cluster>
bool read_clusters(std::istream &stream, std::vector<Cluster> *cluster);
bool write_clusters(std::ostream &stream, const std::vector<Cluster> &cluster);

// MetaMatch::Peak
bool read_peak(std::istream &stream, Peak *peak);
bool write_peak(std::ostream &stream, const Peak &peak);

// std::vector<MetaMatch::Peak>
bool read_peaks(std::istream &stream, std::vector<Peak> *peak);
bool write_peaks(std::ostream &stream, const std::vector<Peak> &peak);

}  // namespace MetaMatch::Serialize

#endif /* METAMATCH_METAMATCHSERIALIZE_HPP */