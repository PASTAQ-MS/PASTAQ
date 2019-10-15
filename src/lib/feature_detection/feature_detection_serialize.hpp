#ifndef FEATUREDETECTION_FEATUREDETECTIONSERIALIZE_HPP
#define FEATUREDETECTION_FEATUREDETECTIONSERIALIZE_HPP

#include <iostream>

#include "feature_detection/feature_detection.hpp"

// This namespace groups the functions used to serialize FeatureDetection data
// structures into a binary stream.
namespace FeatureDetection::Serialize {

// FeatureDetection::Feature
bool read_feature(std::istream &stream, Feature *feature);
bool write_feature(std::ostream &stream, const Feature &feature);

}  // namespace FeatureDetection::Serialize

#endif /* FEATUREDETECTION_FEATUREDETECTIONSERIALIZE_HPP */
