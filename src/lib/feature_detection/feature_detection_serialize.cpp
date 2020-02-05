#include "feature_detection/feature_detection_serialize.hpp"
#include "utils/serialization.hpp"

bool FeatureDetection::Serialize::read_feature(std::istream &stream,
                                               Feature *feature) {
    Serialization::read_uint64(stream, &feature->id);
    Serialization::read_double(stream, &feature->score);
    Serialization::read_double(stream, &feature->average_rt);
    Serialization::read_double(stream, &feature->average_rt_delta);
    Serialization::read_double(stream, &feature->average_rt_sigma);
    Serialization::read_double(stream, &feature->average_mz);
    Serialization::read_double(stream, &feature->average_mz_sigma);
    Serialization::read_double(stream, &feature->total_height);
    Serialization::read_double(stream, &feature->total_volume);
    Serialization::read_double(stream, &feature->max_height);
    Serialization::read_double(stream, &feature->max_volume);
    Serialization::read_double(stream, &feature->monoisotopic_mz);
    Serialization::read_double(stream, &feature->monoisotopic_rt);
    Serialization::read_double(stream, &feature->monoisotopic_height);
    Serialization::read_double(stream, &feature->monoisotopic_volume);
    Serialization::read_int8(stream, &feature->charge_state);
    Serialization::read_vector<uint64_t>(stream, &feature->peak_ids,
                                         Serialization::read_uint64);
    return stream.good();
}

bool FeatureDetection::Serialize::write_feature(std::ostream &stream,
                                                const Feature &feature) {
    Serialization::write_uint64(stream, feature.id);
    Serialization::write_double(stream, feature.score);
    Serialization::write_double(stream, feature.average_rt);
    Serialization::write_double(stream, feature.average_rt_delta);
    Serialization::write_double(stream, feature.average_rt_sigma);
    Serialization::write_double(stream, feature.average_mz);
    Serialization::write_double(stream, feature.average_mz_sigma);
    Serialization::write_double(stream, feature.total_height);
    Serialization::write_double(stream, feature.total_volume);
    Serialization::write_double(stream, feature.max_height);
    Serialization::write_double(stream, feature.max_volume);
    Serialization::write_double(stream, feature.monoisotopic_mz);
    Serialization::write_double(stream, feature.monoisotopic_rt);
    Serialization::write_double(stream, feature.monoisotopic_height);
    Serialization::write_double(stream, feature.monoisotopic_volume);
    Serialization::write_int8(stream, feature.charge_state);
    Serialization::write_vector<uint64_t>(stream, feature.peak_ids,
                                          Serialization::write_uint64);
    return stream.good();
}

bool FeatureDetection::Serialize::read_features(
    std::istream &stream, std::vector<Feature> *features) {
    return Serialization::read_vector<Feature>(
        stream, features, FeatureDetection::Serialize::read_feature);
}

bool FeatureDetection::Serialize::write_features(
    std::ostream &stream, const std::vector<Feature> &features) {
    return Serialization::write_vector<Feature>(
        stream, features, FeatureDetection::Serialize::write_feature);
}
