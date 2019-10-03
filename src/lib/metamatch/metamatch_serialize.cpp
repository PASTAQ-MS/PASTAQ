#include "metamatch/metamatch_serialize.hpp"
#include "centroid/centroid_serialize.hpp"
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

bool MetaMatch::Serialize::read_clusters(std::istream &stream,
                                         std::vector<Cluster> *clusters) {
    return Serialization::read_vector<Cluster>(stream, clusters, read_cluster);
}

bool MetaMatch::Serialize::write_clusters(
    std::ostream &stream, const std::vector<Cluster> &clusters) {
    return Serialization::write_vector<Cluster>(stream, clusters,
                                                write_cluster);
}

bool MetaMatch::Serialize::read_peak(std::istream &stream, Peak *peak) {
    Centroid::Serialize::read_peak(stream, peak);
    Serialization::read_uint32(stream, &peak->file_id);
    Serialization::read_uint32(stream, &peak->class_id);
    Serialization::read_int64(stream, &peak->cluster_id);
    Serialization::read_double(stream, &peak->cluster_mz);
    Serialization::read_double(stream, &peak->cluster_rt);
    return stream.good();
}

bool MetaMatch::Serialize::write_peak(std::ostream &stream, const Peak &peak) {
    Centroid::Serialize::write_peak(stream, peak);
    Serialization::write_uint64(stream, peak.file_id);
    Serialization::write_uint64(stream, peak.class_id);
    Serialization::write_int64(stream, peak.cluster_id);
    Serialization::write_double(stream, peak.cluster_mz);
    Serialization::write_double(stream, peak.cluster_rt);
    return stream.good();
}

bool MetaMatch::Serialize::read_peaks(std::istream &stream,
                                      std::vector<Peak> *peaks) {
    uint64_t num_peaks = 0;
    Serialization::read_uint64(stream, &num_peaks);
    *peaks = std::vector<Peak>(num_peaks);
    for (size_t i = 0; i < num_peaks; ++i) {
        MetaMatch::Serialize::read_peak(stream, &(*peaks)[i]);
    }
    return stream.good();
}

bool MetaMatch::Serialize::write_peaks(std::ostream &stream,
                                       const std::vector<Peak> &peaks) {
    uint64_t num_peaks = peaks.size();
    Serialization::write_uint64(stream, num_peaks);
    for (size_t i = 0; i < num_peaks; ++i) {
        MetaMatch::Serialize::write_peak(stream, peaks[i]);
    }
    return stream.good();
}

