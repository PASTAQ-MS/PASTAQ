#ifndef METAMATCH_METAMATCHSERIALIZE_HPP
#define METAMATCH_METAMATCHSERIALIZE_HPP

#include <iostream>

#include "metamatch.hpp"

// This namespace groups the functions used to serialize MetaMatch data
// structures into a binary stream.
namespace MetaMatch::Serialize {

// MetaMatch::Cluster
bool read_cluster(std::istream &stream, Cluster *cluster);
bool write_cluster(std::ostream &stream, const Cluster &cluster);

}  // namespace MetaMatch::Serialize

#endif /* METAMATCH_METAMATCHSERIALIZE_HPP */
