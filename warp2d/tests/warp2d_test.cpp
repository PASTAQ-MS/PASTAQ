#include <algorithm>
#include <iostream>

#include "doctest.h"
#include "warp2d/warp2d.hpp"

TEST_CASE("DEBUG ALGORITHM") {
    std::vector<Centroid::Peak> target_peaks;
    std::vector<Centroid::Peak> source_peaks;
    auto warped_peaks =
        Warp2D::warp_peaks(target_peaks, source_peaks, 10, 2, 1);
    CHECK(warped_peaks.size() == source_peaks.size());
    CHECK(1 == 0);
}
