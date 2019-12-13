#include "doctest.h"
#include "test_utils.hpp"

#include "metamatch/metamatch.hpp"

TEST_CASE("Clustering of peak lists") {
    // TODO:...
    CHECK(true);
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
        feature.average_rt_sigma += candidate.fitted_sigma_rt;
        feature.average_mz_sigma += candidate.fitted_sigma_mz;
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
    // Three features with small differences.
    std::vector<Centroid::Peak> peaks_a = {
        TestUtils::mock_gaussian_peak(0, 100.0, 400.0, 2000.0, 0.001, 10),
        TestUtils::mock_gaussian_peak(1, 75.0, 401.0, 2000.0, 0.001, 10),
        TestUtils::mock_gaussian_peak(2, 20.0, 402.0, 2000.0, 0.001, 10),
        TestUtils::mock_gaussian_peak(3, 1.0, 403.0, 2000.0, 0.001, 10),
    };
    std::vector<Centroid::Peak> peaks_b = {
        TestUtils::mock_gaussian_peak(0, 110.0, 400.0, 2002.0, 0.001, 10),
        TestUtils::mock_gaussian_peak(1, 85.0, 401.0, 2002.0, 0.001, 10),
        TestUtils::mock_gaussian_peak(2, 30.0, 402.0, 2002.0, 0.001, 10),
        TestUtils::mock_gaussian_peak(3, 11.0, 403.0, 2002.0, 0.001, 10),
    };
    std::vector<Centroid::Peak> peaks_c = {
        TestUtils::mock_gaussian_peak(4, 0.100, 400.0, 2000.0, 0.001, 10),
        TestUtils::mock_gaussian_peak(5, 0.75, 401.0, 2000.0, 0.001, 10),
        TestUtils::mock_gaussian_peak(6, 0.20, 402.0, 2000.0, 0.001, 10),
        TestUtils::mock_gaussian_peak(7, 0.1, 403.0, 2000.0, 0.001, 10),
    };
    std::vector<FeatureDetection::Feature> features_a = {
        mock_feature(0, peaks_a),
        mock_feature(1, peaks_c),
    };
    std::vector<FeatureDetection::Feature> features_b = {
        mock_feature(0, peaks_b),
    };
    std::vector<MetaMatch::InputSetFeatures> input_sets = {
        {0, peaks_a, features_a},
        {0, peaks_b, features_b},
    };
    auto clusters = MetaMatch::find_feature_clusters(input_sets);
    for (const auto &cluster : clusters) {
        CHECK(cluster.id == 0);
        CHECK(TestUtils::compare_double(cluster.mz, 400.678, 3));
        CHECK(TestUtils::compare_double(cluster.rt, 2001));
        CHECK(cluster.feature_ids.size() == 2);
    }
}
