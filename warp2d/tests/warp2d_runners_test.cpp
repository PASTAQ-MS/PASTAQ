#include "doctest.h"

#include "warp2d/warp2d.hpp"
#include "warp2d/warp2d_runners.hpp"

TEST_CASE("Peak warping using warp2d") {
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
    SUBCASE("Parallel execution") {
        {
            auto warped_peaks = Warp2D::Runners::Parallel::run(
                target_peaks, source_peaks, {5, 10, 1000}, 1);
            bool sizes_match = warped_peaks.size() == source_peaks.size();
            CHECK(sizes_match);
            if (sizes_match) {
                for (size_t i = 0; i < warped_peaks.size(); ++i) {
                    CHECK(std::abs(target_peaks[i].rt - warped_peaks[i].rt) <
                          1);
                }
            }
        }
        {
            auto warped_peaks = Warp2D::Runners::Parallel::run(
                target_peaks, source_peaks, {5, 10, 1000}, 4);
            bool sizes_match = warped_peaks.size() == source_peaks.size();
            CHECK(sizes_match);
            if (sizes_match) {
                for (size_t i = 0; i < warped_peaks.size(); ++i) {
                    CHECK(std::abs(target_peaks[i].rt - warped_peaks[i].rt) <
                          1);
                }
            }
        }
        {
            auto warped_peaks = Warp2D::Runners::Parallel::run(
                target_peaks, source_peaks, {5, 10, 1000}, 8);
            bool sizes_match = warped_peaks.size() == source_peaks.size();
            CHECK(sizes_match);
            if (sizes_match) {
                for (size_t i = 0; i < warped_peaks.size(); ++i) {
                    CHECK(std::abs(target_peaks[i].rt - warped_peaks[i].rt) <
                          1);
                }
            }
        }
        {
            auto warped_peaks = Warp2D::Runners::Parallel::run(
                target_peaks, source_peaks, {5, 10, 1000}, 12);
            bool sizes_match = warped_peaks.size() == source_peaks.size();
            CHECK(sizes_match);
            if (sizes_match) {
                for (size_t i = 0; i < warped_peaks.size(); ++i) {
                    CHECK(std::abs(target_peaks[i].rt - warped_peaks[i].rt) <
                          1);
                }
            }
        }
    }
    SUBCASE("Serial execution") {
        auto warped_peaks = Warp2D::Runners::Serial::run(
            target_peaks, source_peaks, {5, 10, 1000});
        bool sizes_match = warped_peaks.size() == source_peaks.size();
        CHECK(sizes_match);
        if (sizes_match) {
            for (size_t i = 0; i < warped_peaks.size(); ++i) {
                CHECK(std::abs(target_peaks[i].rt - warped_peaks[i].rt) < 1);
            }
        }
    }
}
