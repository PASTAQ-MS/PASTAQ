#include "doctest.h"
#include "test_utils.hpp"

#include "warp2d/warp2d.hpp"

TEST_CASE("Peak warping using warp2d") {
    std::vector<Centroid::Peak> reference_peaks = {};
    std::vector<Centroid::Peak> peaks = {};
    SUBCASE("Parallel execution (1 thread)") {
        // TODO:...
        CHECK(true);
    }
    SUBCASE("Parallel execution (2 threads)") {
        // TODO:...
        CHECK(true);
    }
    SUBCASE("Parallel execution (4 threads)") {
        // TODO:...
        CHECK(true);
    }
    SUBCASE("Serial execution") {
        // TODO:...
        CHECK(true);
    }
}
