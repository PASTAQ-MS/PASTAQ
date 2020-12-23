#include "metamatch/metamatch_serialize.hpp"
#include "utils/serialization.hpp"

bool read_peak_id(std::istream &stream, MetaMatch::PeakId *peak_id) {
    Serialization::read_uint64(stream, &peak_id->file_id);
    Serialization::read_uint64(stream, &peak_id->peak_id);
    return stream.good();
}

bool write_peak_id(std::ostream &stream, const MetaMatch::PeakId &peak_id) {
    Serialization::write_uint64(stream, peak_id.file_id);
    Serialization::write_uint64(stream, peak_id.peak_id);
    return stream.good();
}

bool read_peak_cluster(std::istream &stream, MetaMatch::PeakCluster *cluster) {
    Serialization::read_uint64(stream, &cluster->id);
    Serialization::read_double(stream, &cluster->mz);
    Serialization::read_double(stream, &cluster->rt);
    Serialization::read_double(stream, &cluster->avg_height);
    Serialization::read_double(stream, &cluster->avg_volume);

    Serialization::read_vector<double>(stream, &cluster->heights,
                                       Serialization::read_double);
    Serialization::read_vector<double>(stream, &cluster->volumes,
                                       Serialization::read_double);

    Serialization::read_vector<MetaMatch::PeakId>(stream, &cluster->peak_ids,
                                          read_peak_id);
    return stream.good();
}

bool write_peak_cluster(
    std::ostream &stream,  const MetaMatch::PeakCluster &cluster) {
    Serialization::write_uint64(stream, cluster.id);
    Serialization::write_double(stream, cluster.mz);
    Serialization::write_double(stream, cluster.rt);
    Serialization::write_double(stream, cluster.avg_height);
    Serialization::write_double(stream, cluster.avg_volume);

    Serialization::write_vector<double>(stream, cluster.heights,
                                        Serialization::write_double);
    Serialization::write_vector<double>(stream, cluster.volumes,
                                        Serialization::write_double);

    Serialization::write_vector<MetaMatch::PeakId>(stream, cluster.peak_ids, write_peak_id);
    return stream.good();
}

bool MetaMatch::Serialize::read_peak_clusters(
    std::istream &stream, std::vector<MetaMatch::PeakCluster> *clusters) {
    Serialization::read_vector<MetaMatch::PeakCluster>(stream, clusters,
                                               read_peak_cluster);
    return stream.good();
}

bool MetaMatch::Serialize::write_peak_clusters(
    std::ostream &stream, const std::vector<PeakCluster> &clusters) {
    Serialization::write_vector<PeakCluster>(stream, clusters, write_peak_cluster);
    return stream.good();
}

bool read_feature_id(std::istream &stream, MetaMatch::FeatureId *feature_id) {
    Serialization::read_uint64(stream, &feature_id->file_id);
    Serialization::read_uint64(stream, &feature_id->feature_id);
    return stream.good();
}

bool write_feature_id(std::ostream &stream, const MetaMatch::FeatureId &feature_id) {
    Serialization::write_uint64(stream, feature_id.file_id);
    Serialization::write_uint64(stream, feature_id.feature_id);
    return stream.good();
}

bool read_feature_cluster(std::istream &stream, MetaMatch::FeatureCluster *cluster) {
    Serialization::read_uint64(stream, &cluster->id);
    Serialization::read_double(stream, &cluster->mz);
    Serialization::read_double(stream, &cluster->rt);
    Serialization::read_int8(stream, &cluster->charge_state);
    Serialization::read_double(stream, &cluster->avg_total_height);
    Serialization::read_double(stream, &cluster->avg_monoisotopic_height);
    Serialization::read_double(stream, &cluster->avg_max_height);
    Serialization::read_double(stream, &cluster->avg_total_volume);
    Serialization::read_double(stream, &cluster->avg_monoisotopic_volume);
    Serialization::read_double(stream, &cluster->avg_max_volume);

    Serialization::read_vector<double>(stream, &cluster->total_heights,
                                       Serialization::read_double);
    Serialization::read_vector<double>(stream, &cluster->monoisotopic_heights,
                                       Serialization::read_double);
    Serialization::read_vector<double>(stream, &cluster->max_heights,
                                       Serialization::read_double);
    Serialization::read_vector<double>(stream, &cluster->total_volumes,
                                       Serialization::read_double);
    Serialization::read_vector<double>(stream, &cluster->monoisotopic_volumes,
                                       Serialization::read_double);
    Serialization::read_vector<double>(stream, &cluster->max_volumes,
                                       Serialization::read_double);

    Serialization::read_vector<MetaMatch::FeatureId>(stream, &cluster->feature_ids,
                                          read_feature_id);
    return stream.good();
}

bool write_feature_cluster(
    std::ostream &stream,  const MetaMatch::FeatureCluster &cluster) {
    Serialization::write_uint64(stream, cluster.id);
    Serialization::write_double(stream, cluster.mz);
    Serialization::write_double(stream, cluster.rt);
    Serialization::write_int8(stream, cluster.charge_state);
    Serialization::write_double(stream, cluster.avg_total_height);
    Serialization::write_double(stream, cluster.avg_monoisotopic_height);
    Serialization::write_double(stream, cluster.avg_max_height);
    Serialization::write_double(stream, cluster.avg_total_volume);
    Serialization::write_double(stream, cluster.avg_monoisotopic_volume);
    Serialization::write_double(stream, cluster.avg_max_volume);

    Serialization::write_vector<double>(stream, cluster.total_heights,
                                        Serialization::write_double);
    Serialization::write_vector<double>(stream, cluster.monoisotopic_heights,
                                        Serialization::write_double);
    Serialization::write_vector<double>(stream, cluster.max_heights,
                                        Serialization::write_double);
    Serialization::write_vector<double>(stream, cluster.total_volumes,
                                        Serialization::write_double);
    Serialization::write_vector<double>(stream, cluster.monoisotopic_volumes,
                                        Serialization::write_double);
    Serialization::write_vector<double>(stream, cluster.max_volumes,
                                        Serialization::write_double);

    Serialization::write_vector<MetaMatch::FeatureId>(stream, cluster.feature_ids, write_feature_id);
    return stream.good();
}

bool MetaMatch::Serialize::read_feature_clusters(
    std::istream &stream, std::vector<MetaMatch::FeatureCluster> *clusters) {
    Serialization::read_vector<MetaMatch::FeatureCluster>(stream, clusters,
                                               read_feature_cluster);
    return stream.good();
}

bool MetaMatch::Serialize::write_feature_clusters(
    std::ostream &stream, const std::vector<FeatureCluster> &clusters) {
    Serialization::write_vector<FeatureCluster>(stream, clusters, write_feature_cluster);
    return stream.good();
}
