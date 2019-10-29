#ifndef FEATUREDETECTION__FEATUREDETECTION_HPP
#define FEATUREDETECTION__FEATUREDETECTION_HPP
#include <optional>

#include "link/link.hpp"
#include "utils/search.hpp"

namespace FeatureDetection {
// TODO: Expand documentation.
struct Feature {
    uint64_t id;
    int64_t msms_id;  // (Optional, might not exist if value is < 0)
                      // TODO: Should we make this an std::optional type?
                      // TODO: Should we create a 'link table' for features to
                      // not pollute the feature namespace?
    double average_rt;
    double average_rt_delta;
    double average_rt_sigma;
    double average_mz;
    double average_mz_sigma;
    double total_height;
    double monoisotopic_mz;
    double monoisotopic_height;
    // FIXME: Use the peak index instead of the peak_id? Tradeoffs?
    std::vector<uint64_t> peak_ids;
};

// TODO: Expand documentation.
struct TheoreticalIsotopes {
    std::vector<double> mzs;
    std::vector<double> percs;
};

TheoreticalIsotopes theoretical_isotopes_peptide(std::string sequence,
                                                 int8_t charge_state,
                                                 double min_perc);

// TODO: Too many parameters, probably we want to pass a parameters struct
// instead.
std::optional<Feature> build_feature(
    const std::vector<bool> &peaks_in_use,
    const std::vector<Centroid::Peak> &peaks,
    const std::vector<Search::KeySort<double>> &peaks_rt_key,
    const TheoreticalIsotopes &theoretical_isotopes, double tolerance_rt,
    double retention_time, double discrepancy_threshold);

std::vector<Feature> feature_detection(
    const std::vector<Centroid::Peak> &peaks,
    const RawData::RawData &raw_data_ms2,
    const IdentData::IdentData &ident_data,
    const std::vector<Link::LinkedMsms> &link_table_msms,
    const std::vector<Link::LinkedMsms> &link_table_idents,
    double discrepancy_threshold);

}  // namespace FeatureDetection

#endif /* FEATUREDETECTION__FEATUREDETECTION_HPP */
