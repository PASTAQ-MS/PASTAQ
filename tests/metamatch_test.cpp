#include "doctest.h"
#include "test_utils.hpp"

#include "metamatch/metamatch.hpp"

TEST_CASE("Clustering of peak lists") {
    // TODO:...
    CHECK(true);
}

Centroid::Peak mock_gaussian_peak(size_t id, double height, double mz,
                                  double rt, double sigma_mz, double sigma_rt) {
    Centroid::Peak peak = {};
    peak.id = id;
    peak.local_max_mz = mz;
    peak.local_max_rt = rt;
    peak.local_max_height = height;
    peak.rt_delta = 0.0;
    peak.roi_min_mz = mz - 3 * sigma_mz;
    peak.roi_max_mz = mz + 3 * sigma_mz;
    peak.roi_min_rt = rt - 3 * sigma_rt;
    peak.roi_max_rt = rt + 3 * sigma_rt;
    peak.raw_roi_mean_mz = mz;
    peak.raw_roi_mean_rt = rt;
    peak.raw_roi_sigma_mz = sigma_mz;
    peak.raw_roi_sigma_rt = sigma_rt;
    peak.raw_roi_max_height = height;
    peak.raw_roi_total_intensity = height;
    return peak;
}

FeatureDetection::Feature mock_feature(size_t id,
                                       std::vector<Centroid::Peak> &peaks) {
    FeatureDetection::Feature feature = {};
    feature.id = id;
    feature.monoisotopic_mz = peaks[0].local_max_mz;
    feature.monoisotopic_height = peaks[0].local_max_height;
    // Find the weighted average mz and the average height of the selected
    // isotopes.
    feature.average_mz = 0.0;
    feature.total_height = 0.0;
    feature.average_rt = 0.0;
    feature.average_rt_delta = 0.0;
    feature.average_rt_sigma = 0.0;
    feature.average_mz_sigma = 0.0;
    for (size_t i = 0; i < peaks.size(); ++i) {
        auto candidate = peaks[i];
        feature.total_height += candidate.local_max_height;
        feature.average_mz +=
            candidate.local_max_height * candidate.local_max_mz;
        feature.average_rt += candidate.local_max_rt;
        feature.average_rt_delta += candidate.rt_delta;
        feature.average_rt_sigma += candidate.raw_roi_sigma_rt;
        feature.average_mz_sigma += candidate.raw_roi_sigma_mz;
        feature.peak_ids.push_back(candidate.id);
    }
    feature.average_mz /= feature.total_height;
    feature.average_rt /= peaks.size();
    feature.average_rt_delta /= peaks.size();
    feature.average_rt_sigma /= peaks.size();
    feature.average_mz_sigma /= peaks.size();
    return feature;
}

TEST_CASE("Clustering features lists") {
    // TODO:...
    // Three features with small differences.
    std::vector<Centroid::Peak> peaks_a = {
        mock_gaussian_peak(0, 100.0, 400.0, 2000.0, 0.001, 10),
        mock_gaussian_peak(1, 75.0, 401.0, 2000.0, 0.001, 10),
        mock_gaussian_peak(2, 20.0, 402.0, 2000.0, 0.001, 10),
        mock_gaussian_peak(3, 1.0, 403.0, 2000.0, 0.001, 10),
    };
    std::vector<Centroid::Peak> peaks_b = {
        mock_gaussian_peak(0, 100.0, 400.0, 2002.0, 0.001, 10),
        mock_gaussian_peak(1, 75.0, 401.0, 2002.0, 0.001, 10),
        mock_gaussian_peak(2, 20.0, 402.0, 2002.0, 0.001, 10),
        mock_gaussian_peak(3, 1.0, 403.0, 2002.0, 0.001, 10),
    };
    std::vector<FeatureDetection::Feature> features_a = {
        mock_feature(0, peaks_a),
    };
    std::vector<FeatureDetection::Feature> features_b = {
        mock_feature(0, peaks_b),
    };
    std::vector<MetaMatch::InputSetFeatures> input_sets = {
        {0, &peaks_a, &features_a},
        {0, &peaks_b, &features_b},
    };
    auto clusters = MetaMatch::find_feature_clusters(input_sets);
    // CHECK(true);
    CHECK(false);
}
