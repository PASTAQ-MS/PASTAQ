#ifndef FEATUREDETECTION_FEATUREDETECTIONSERIALIZE_HPP
#define FEATUREDETECTION_FEATUREDETECTIONSERIALIZE_HPP

#include <iostream>

#include "feature_detection/feature_detection.hpp"

// This namespace groups the functions used to serialize FeatureDetection data
// structures into a binary stream.
namespace FeatureDetection::Serialize {

// Read/write a feature to the given binary stream.
bool read_feature(std::istream &stream, Feature *feature);
bool write_feature(std::ostream &stream, const Feature &feature);

// Read/write all features to the given binary stream.
bool read_features(std::istream &stream, std::vector<Feature> *features);
bool write_features(std::ostream &stream, const std::vector<Feature> &features);

}  // namespace FeatureDetection::Serialize

#endif /* FEATUREDETECTION_FEATUREDETECTIONSERIALIZE_HPP */
