#include <algorithm>
#include <iostream>

#include "doctest.h"
#include "warp2d/warp2d.hpp"

TEST_CASE("DEBUG ALGORITHM") {
    std::vector<Centroid::Peak> target_peaks = {
        {0, 0, 100, 300, 1000, 10000, 5, 200, 100, 300, 1000, 10000, 0, {}, {}},
        {0, 0, 120, 300, 900, 9000, 5, 200, 120, 300, 900, 9000, 0, {}, {}},
        {0, 0, 140, 300, 700, 7000, 5, 200, 140, 300, 700, 7000, 0, {}, {}},
        {0, 0, 200, 600, 1000, 10000, 5, 200, 200, 600, 1000, 10000, 0, {}, {}},
        {0, 0, 220, 600, 900, 9000, 5, 200, 220, 600, 900, 9000, 0, {}, {}},
        {0, 0, 240, 600, 700, 7000, 5, 200, 240, 600, 700, 7000, 0, {}, {}},
    };
    // Previous peaks shifted 50 units down.
    std::vector<Centroid::Peak> source_peaks = {
        {0, 0, 100, 350, 1000, 10000, 5, 200, 100, 350, 1000, 10000, 0, {}, {}},
        {0, 0, 120, 350, 900, 9000, 5, 200, 120, 350, 900, 9000, 0, {}, {}},
        {0, 0, 140, 350, 700, 7000, 5, 200, 140, 350, 700, 7000, 0, {}, {}},
        {0, 0, 200, 650, 1000, 10000, 5, 200, 200, 650, 1000, 10000, 0, {}, {}},
        {0, 0, 220, 650, 900, 9000, 5, 200, 220, 650, 900, 9000, 0, {}, {}},
        {0, 0, 240, 650, 700, 7000, 5, 200, 240, 650, 700, 7000, 0, {}, {}},
    };
    auto warped_peaks =
        Warp2D::warp_peaks(target_peaks, source_peaks, 10, 2, 1);
    CHECK(warped_peaks.size() == source_peaks.size());
    CHECK(1 == 0);
}

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
