#include "doctest.h"

#include "warp2d/warp2d.hpp"

TEST_CASE("Peak overlapping") {
    auto round_double = [](double d) {
        return (long long int)(d * 1000.0) / 1000.0;
    };

    {
        Centroid::Peak peak_a = {0,   0,   100,  300,   1000, 10000, 5, 200,
                                 100, 300, 1000, 10000, 0,    {},    {}};
        Centroid::Peak peak_b = {0,   0,   100,  300,   1000, 10000, 5, 200,
                                 100, 300, 1000, 10000, 0,    {},    {}};
        double overlap = Warp2D::peak_overlap(peak_a, peak_b);
        CHECK(round_double(500.0) == round_double(overlap));
    }
    {
        Centroid::Peak peak_a = {0,   0,   100,  300,   1000, 10000, 5, 200,
                                 100, 300, 1000, 10000, 0,    {},    {}};
        Centroid::Peak peak_b = {0,   0,   100,  350,   1000, 10000, 5, 200,
                                 100, 350, 1000, 10000, 0,    {},    {}};
        double overlap = Warp2D::peak_overlap(peak_a, peak_b);
        CHECK(round_double(492.248) == round_double(overlap));
    }
    {
        Centroid::Peak peak_a = {0,   0,   100,  300,   1000, 10000, 5, 200,
                                 100, 300, 1000, 10000, 0,    {},    {}};
        Centroid::Peak peak_b = {0,   0,   120, 300,  900, 9000, 5, 200,
                                 120, 300, 900, 9000, 0,   {},   {}};
        double overlap = Warp2D::peak_overlap(peak_a, peak_b);
        CHECK(round_double(8.24204) == round_double(overlap));
    }
    {
        Centroid::Peak peak_a = {0,   0,   100,  300,   1000, 10000, 5, 200,
                                 100, 300, 1000, 10000, 0,    {},    {}};
        Centroid::Peak peak_b = {0,   0,   240, 650,  700, 7000, 5, 200,
                                 240, 650, 700, 7000, 0,   {},   {}};
        double overlap = Warp2D::peak_overlap(peak_a, peak_b);
        CHECK(round_double(0) == round_double(overlap));
    }
}
