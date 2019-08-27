#include "metamatch_serialize.hpp"
#include "utils/serialization.hpp"

bool MetaMatch::Serialize::read_cluster(std::istream &stream,
                                        Cluster *cluster) {
    Serialization::read_int64(stream, &cluster->id);
    Serialization::read_double(stream, &cluster->mz);
    Serialization::read_double(stream, &cluster->rt);
    Serialization::read_double(stream, &cluster->avg_height);
    uint64_t num_files = 0;
    Serialization::read_uint64(stream, &num_files);
    cluster->file_heights = std::vector<double>(num_files);
    for (size_t i = 0; i < num_files; ++i) {
        Serialization::read_double(stream, &cluster->file_heights[i]);
    }
    return stream.good();
}

bool MetaMatch::Serialize::write_cluster(std::ostream &stream,
                                         const Cluster &cluster) {
    Serialization::write_int64(stream, cluster.id);
    Serialization::write_double(stream, cluster.mz);
    Serialization::write_double(stream, cluster.rt);
    Serialization::write_double(stream, cluster.avg_height);
    uint64_t num_files = cluster.file_heights.size();
    Serialization::write_uint64(stream, num_files);
    for (size_t i = 0; i < num_files; ++i) {
        Serialization::write_double(stream, cluster.file_heights[i]);
    }
    return stream.good();
}
