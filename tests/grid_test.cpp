#include "doctest.h"
#include "test_utils.hpp"

#include "grid/grid.hpp"

TEST_CASE("Gaussian splatting") {
    SUBCASE("Splat on the center of the mesh") {
        // TODO:...
        CHECK(true);
    }

    SUBCASE("Splatting outside the grid") {
        // TODO:...
        CHECK(true);
    }

    SUBCASE("Splat outside but range falls inside") {
        // TODO:...
        CHECK(true);
    }
}

TEST_CASE("Test the warped grid interpolation") {
    SUBCASE("QUAD") {
        // TODO:...
        CHECK(true);
    }
    SUBCASE("ORBITRAP") {
        // TODO:...
        CHECK(true);
    }
    SUBCASE("FTICR") {
        // TODO:...
        CHECK(true);
    }
    SUBCASE("TOF") {
        // TODO:...
        CHECK(true);
    }
}

TEST_CASE("Parallel and serial execution offer the same results") {
    // TODO:...
    CHECK(true);
}
