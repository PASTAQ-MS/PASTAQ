#ifndef METAMATCH_METAMATCHFILES_HPP
#define METAMATCH_METAMATCHFILES_HPP
#include <iostream>

#include "metamatch/metamatch.hpp"

namespace MetaMatch::Files::Csv {
// Write the most significant peak values on the csv format.
bool write_clusters(std::ostream& stream,
                    const std::vector<MetaMatch::Cluster>& clusters,
                    size_t n_files);
// Write the peaks to the given stream.
bool write_peaks(std::ostream& stream,
                 const std::vector<MetaMatch::Peak>& peaks,
                 bool include_mpid=false);
}  // namespace MetaMatch::Files::Csv

#endif /* METAMATCH_METAMATCHFILES_HPP */
