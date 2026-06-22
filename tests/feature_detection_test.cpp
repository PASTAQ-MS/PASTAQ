#include "doctest.h"
#include "test_utils.hpp"

#include "feature_detection/feature_detection.hpp"

TEST_CASE("feature detection smoke") {
    std::vector<Centroid::Peak> peaks = {
        TestUtils::mock_gaussian_peak(0, 100.0, 400.0, 2000.0, 0.01, 10),
        TestUtils::mock_gaussian_peak(1, 45.0, 400.5, 2000.0, 0.01, 10),
        TestUtils::mock_gaussian_peak(2, 12.0, 401.0, 2000.0, 0.01, 10),
        TestUtils::mock_gaussian_peak(3, 2.0, 401.5, 2000.0, 0.01, 10),
        TestUtils::mock_gaussian_peak(4, 100.0, 400.5, 2000.0, 0.01, 10),
        TestUtils::mock_gaussian_peak(5, 22.0, 401.5, 2000.0, 0.01, 10),
        TestUtils::mock_gaussian_peak(6, 3.0, 402.5, 2000.0, 0.01, 10),
        TestUtils::mock_gaussian_peak(7, 300.0, 402.0, 2000.0, 0.01, 10),
        TestUtils::mock_gaussian_peak(8, 100.0, 402.0, 2000.0, 0.01, 10),
        TestUtils::mock_gaussian_peak(9, 45.0, 402.5, 2000.0, 0.01, 10),
        TestUtils::mock_gaussian_peak(10, 12.0, 403.0, 2000.0, 0.01, 10),
        TestUtils::mock_gaussian_peak(11, 2.0, 403.5, 2000.0, 0.01, 10),
        TestUtils::mock_gaussian_peak(12, 30000.0, 402.0, 2000.0, 0.01, 10),
        TestUtils::mock_gaussian_peak(13, 30000.0, 402.5, 2000.0, 0.01, 10),
    };
    std::vector<uint8_t> charge_states = {2, 1};

    CHECK_NOTHROW(FeatureDetection::detect_features(peaks, charge_states));
}

TEST_CASE("rolling cosine selects expected window and score") {
    std::vector<double> path_heights = {1.0, 100.0, 50.0, 25.0, 1.0};
    std::vector<double> averagine_heights = {100.0, 50.0, 25.0};

    auto sim = FeatureDetection::Detail::rolling_cosine_sim(path_heights,
                                                             averagine_heights);

    CHECK(sim.min_i == 1);
    CHECK(sim.max_i == 3);
    CHECK(sim.dot == doctest::Approx(1.0).epsilon(1e-12));
}
