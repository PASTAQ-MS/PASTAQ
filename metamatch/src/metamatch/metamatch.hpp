#ifndef METAMATCH_METAMATCH_HPP
#define METAMATCH_METAMATCH_HPP

#include <cstdint>
#include <vector>

#include "centroid/centroid_files.hpp"

namespace MetaMatch {

// TODO(alex): Add docs...
struct Parameters {
    double radius_mz;
    double radius_rt;
    double fraction;
};

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
void find_candidates(std::vector<MetaMatch::Peak>& peaks,
                     const MetaMatch::Parameters& parameters);

// TODO(alex): Add docs...
std::vector<MetaMatch::Peak> extract_orphans(
    std::vector<MetaMatch::Peak>& peaks);

// TODO(alex): Add docs...
std::vector<MetaMatch::Cluster> reduce_cluster(
    std::vector<MetaMatch::Peak>& peaks, size_t n_files);

// TODO(alex): Move to MetaMatch::Files::Csv::write_clusters
bool write_clusters(std::ostream& stream,
                    const std::vector<MetaMatch::Cluster>& clusters,
                    size_t n_files);
// TODO(alex): Move to MetaMatch::Files::Csv::read_file_list
std::vector<std::pair<std::string, size_t>> read_file_list(
    std::istream& stream);
}  // namespace MetaMatch

#endif /* METAMATCH_METAMATCH_HPP */
