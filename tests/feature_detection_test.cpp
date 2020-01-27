#include "doctest.h"
#include "test_utils.hpp"

#include "feature_detection/feature_detection.hpp"

TEST_CASE("DEBUG") {
    std::vector<Centroid::Peak> peaks = {
        TestUtils::mock_gaussian_peak(0, 100.0, 400.0, 2000.0, 0.001, 10),
        TestUtils::mock_gaussian_peak(1, 45.0, 400.5, 2000.0, 0.001, 10),
        TestUtils::mock_gaussian_peak(2, 12.0, 401.0, 2000.0, 0.001, 10),
        TestUtils::mock_gaussian_peak(3, 2.0, 401.5, 2000.0, 0.001, 10),
        TestUtils::mock_gaussian_peak(4, 100.0, 400.5, 2000.0, 0.001, 10),
        TestUtils::mock_gaussian_peak(5, 22.0, 401.5, 2000.0, 0.001, 10),
        TestUtils::mock_gaussian_peak(6, 3.0, 402.5, 2000.0, 0.001, 10),
        TestUtils::mock_gaussian_peak(7, 101010.0, 402.0, 2000.0, 0.001, 10),
        TestUtils::mock_gaussian_peak(8, 101001.0, 402.5, 2000.0, 0.001, 10),
    };
    std::vector<uint8_t> charge_states = {2, 1};
    FeatureDetection::find_candidates(peaks, charge_states);
    CHECK(false);
}

