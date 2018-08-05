#include <algorithm>
#include <iostream>

#include "doctest.h"
#include "warp2d/warp2d.hpp"

void print_peaks_debug(std::vector<Centroid::Peak> &peaks) {
    for (const auto &peak : peaks) {
        std::cout << "rt: " << peak.rt << std::endl;
        std::cout << "mz: " << peak.mz << std::endl;
        std::cout << "height: " << peak.height << std::endl;
        std::cout << "rt_centroid: " << peak.rt_centroid << std::endl;
        std::cout << "mz_centroid: " << peak.mz_centroid << std::endl;
        std::cout << "height_centroid: " << peak.height_centroid << std::endl;
    }
}

TEST_CASE("DEBUG") {
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
        {0, 0, 100, 305, 1000, 10000, 5, 200, 100, 305, 1000, 10000, 0, {}, {}},
        {0, 0, 120, 305, 900, 9000, 5, 200, 120, 305, 900, 9000, 0, {}, {}},
        {0, 0, 140, 305, 700, 7000, 5, 200, 140, 305, 700, 7000, 0, {}, {}},
        {0, 0, 200, 610, 1000, 10000, 5, 200, 200, 610, 1000, 10000, 0, {}, {}},
        {0, 0, 220, 610, 900, 9000, 5, 200, 220, 610, 900, 9000, 0, {}, {}},
        {0, 0, 240, 610, 700, 7000, 5, 200, 240, 610, 700, 7000, 0, {}, {}},
    };
    auto warped_peaks =
        Warp2D::warp_peaks(target_peaks, source_peaks, {1, 5, 25});
    bool sizes_match = warped_peaks.size() == source_peaks.size();
    CHECK(sizes_match);
    //if (sizes_match) {
        //for (size_t i = 0; i < warped_peaks.size(); ++i) {
            //CHECK(std::abs(target_peaks[i].rt - warped_peaks[i].rt) < 1);
        //}
    //}
    CHECK(false);
}

//TEST_CASE("Peak warping using warp2d") {
    //std::vector<Centroid::Peak> target_peaks = {
        //{0, 0, 100, 300, 1000, 10000, 5, 200, 100, 300, 1000, 10000, 0, {}, {}},
        //{0, 0, 120, 300, 900, 9000, 5, 200, 120, 300, 900, 9000, 0, {}, {}},
        //{0, 0, 140, 300, 700, 7000, 5, 200, 140, 300, 700, 7000, 0, {}, {}},
        //{0, 0, 200, 600, 1000, 10000, 5, 200, 200, 600, 1000, 10000, 0, {}, {}},
        //{0, 0, 220, 600, 900, 9000, 5, 200, 220, 600, 900, 9000, 0, {}, {}},
        //{0, 0, 240, 600, 700, 7000, 5, 200, 240, 600, 700, 7000, 0, {}, {}},
    //};
    //// Previous peaks shifted 50 units down.
    //std::vector<Centroid::Peak> source_peaks = {
        //{0, 0, 100, 305, 1000, 10000, 5, 200, 100, 305, 1000, 10000, 0, {}, {}},
        //{0, 0, 120, 305, 900, 9000, 5, 200, 120, 305, 900, 9000, 0, {}, {}},
        //{0, 0, 140, 305, 700, 7000, 5, 200, 140, 305, 700, 7000, 0, {}, {}},
        //{0, 0, 200, 610, 1000, 10000, 5, 200, 200, 610, 1000, 10000, 0, {}, {}},
        //{0, 0, 220, 610, 900, 9000, 5, 200, 220, 610, 900, 9000, 0, {}, {}},
        //{0, 0, 240, 610, 700, 7000, 5, 200, 240, 610, 700, 7000, 0, {}, {}},
    //};
    //auto warped_peaks =
        //Warp2D::warp_peaks(target_peaks, source_peaks, {5, 10, 100});
    //bool sizes_match = warped_peaks.size() == source_peaks.size();
    //CHECK(sizes_match);
    //if (sizes_match) {
        //for (size_t i = 0; i < warped_peaks.size(); ++i) {
            //CHECK(std::abs(target_peaks[i].rt - warped_peaks[i].rt) < 1);
        //}
    //}
//}

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
