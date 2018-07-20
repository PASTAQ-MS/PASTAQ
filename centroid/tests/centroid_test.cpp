#include <algorithm>
#include <iostream>

#include "centroid/centroid.hpp"
#include "doctest.h"

TEST_CASE("Find local maxima") {
    Grid::Parameters grid_params = {
        {}, {0, 100, 0, 100}, {50.0, 5.0, 5.0}, Instrument::QUAD, 0x00};
    Grid::calculate_dimensions(grid_params);
    std::vector<double> data(grid_params.dimensions.n *
                             grid_params.dimensions.m);
    // Simulating one isotope with three adjacent peaks and an artifact peak.
    // Peak 1
    CHECK(Grid::splat({20, 40, 3}, grid_params, data));
    CHECK(Grid::splat({20, 50, 5}, grid_params, data));
    CHECK(Grid::splat({20, 57, 5.5}, grid_params, data));
    CHECK(Grid::splat({20, 60, 6}, grid_params, data));
    CHECK(Grid::splat({20, 63, 5}, grid_params, data));
    CHECK(Grid::splat({20, 65, 4}, grid_params, data));
    // Peak 2
    CHECK(Grid::splat({50, 50, 4}, grid_params, data));
    CHECK(Grid::splat({50, 57, 4.5}, grid_params, data));
    CHECK(Grid::splat({50, 60, 6}, grid_params, data));
    CHECK(Grid::splat({50, 63, 4}, grid_params, data));
    CHECK(Grid::splat({50, 65, 3}, grid_params, data));
    // Peak 3
    CHECK(Grid::splat({80, 57, 4.5}, grid_params, data));
    CHECK(Grid::splat({80, 60, 5}, grid_params, data));
    CHECK(Grid::splat({80, 63, 3}, grid_params, data));
    CHECK(Grid::splat({80, 65, 2}, grid_params, data));
    // Artifact peak
    CHECK(Grid::splat({10, 10, 4}, grid_params, data));

    auto points = Centroid::find_local_maxima({0, 0, grid_params}, data);
    CHECK(points.size() == 4);
}

TEST_CASE("Find peak boundaries") {
    Grid::Parameters grid_params = {
        {}, {0, 100, 0, 100}, {50.0, 5.0, 5.0}, Instrument::QUAD, 0x00};
    Grid::calculate_dimensions(grid_params);
    std::vector<double> data(grid_params.dimensions.n *
                             grid_params.dimensions.m);
    // Simulating one isotope with three adjacent peaks and an artifact peak.
    // Peak 1
    CHECK(Grid::splat({20, 40, 3}, grid_params, data));
    CHECK(Grid::splat({20, 50, 5}, grid_params, data));
    CHECK(Grid::splat({20, 57, 5.5}, grid_params, data));
    CHECK(Grid::splat({20, 60, 6}, grid_params, data));
    CHECK(Grid::splat({20, 63, 5}, grid_params, data));
    CHECK(Grid::splat({20, 65, 4}, grid_params, data));
    // Peak 2
    CHECK(Grid::splat({50, 50, 4}, grid_params, data));
    CHECK(Grid::splat({50, 57, 4.5}, grid_params, data));
    CHECK(Grid::splat({50, 60, 6}, grid_params, data));
    CHECK(Grid::splat({50, 63, 4}, grid_params, data));
    CHECK(Grid::splat({50, 65, 3}, grid_params, data));
    // Peak 3
    CHECK(Grid::splat({80, 57, 4.5}, grid_params, data));
    CHECK(Grid::splat({80, 60, 5}, grid_params, data));
    CHECK(Grid::splat({80, 63, 3}, grid_params, data));
    CHECK(Grid::splat({80, 65, 2}, grid_params, data));
    // Artifact peak
    CHECK(Grid::splat({10, 10, 4}, grid_params, data));

    auto local_max_points =
        Centroid::find_local_maxima({0, 0, grid_params}, data);
    CHECK(local_max_points.size() == 4);

    std::vector<size_t> expected_size = {
        50,
        40,
        35,
        25,
    };

    // Only interested in the indices, ignoring the value for now, but we assume
    // it's correct.
    std::vector<std::vector<Centroid::Point>> expected_boundary = {
        {
            {2, 6, 0},  {3, 6, 0},  {4, 6, 0},  {5, 6, 0},  {6, 6, 0},
            {2, 7, 0},  {6, 7, 0},  {2, 8, 0},  {6, 8, 0},  {2, 9, 0},
            {6, 9, 0},  {2, 10, 0}, {6, 10, 0}, {2, 11, 0}, {6, 11, 0},
            {2, 12, 0}, {6, 12, 0}, {2, 13, 0}, {6, 13, 0}, {2, 14, 0},
            {6, 14, 0}, {2, 15, 0}, {3, 15, 0}, {4, 15, 0}, {5, 15, 0},
            {6, 15, 0},
        },
        {
            {8, 8, 0},   {9, 8, 0},   {10, 8, 0},  {11, 8, 0},  {12, 8, 0},
            {8, 9, 0},   {12, 9, 0},  {8, 10, 0},  {12, 10, 0}, {8, 11, 0},
            {12, 11, 0}, {8, 12, 0},  {12, 12, 0}, {8, 13, 0},  {12, 13, 0},
            {8, 14, 0},  {12, 14, 0}, {8, 15, 0},  {9, 15, 0},  {10, 15, 0},
            {11, 15, 0}, {12, 15, 0},
        },
        {
            {14, 9, 0},  {15, 9, 0},  {16, 9, 0},  {17, 9, 0},  {18, 9, 0},
            {14, 10, 0}, {18, 10, 0}, {14, 11, 0}, {18, 11, 0}, {14, 12, 0},
            {18, 12, 0}, {14, 13, 0}, {18, 13, 0}, {14, 14, 0}, {18, 14, 0},
            {14, 15, 0}, {15, 15, 0}, {16, 15, 0}, {17, 15, 0}, {18, 15, 0},
        },
        {
            {0, 0, 0},
            {1, 0, 0},
            {2, 0, 0},
            {3, 0, 0},
            {4, 0, 0},
            {0, 1, 0},
            {4, 1, 0},
            {0, 2, 0},
            {4, 2, 0},
            {0, 3, 0},
            {4, 3, 0},
            {0, 4, 0},
            {1, 4, 0},
            {2, 4, 0},
            {3, 4, 0},
            {4, 4, 0},
        },
    };
    for (size_t i = 0; i < local_max_points.size(); ++i) {
        auto local_max = local_max_points[i];
        std::vector<Centroid::Point> peak_points;
        Centroid::explore_peak_slope(local_max.i, local_max.j, -1,
                                     {0, 0, grid_params}, data, peak_points);
        CHECK(expected_size[i] == peak_points.size());
        auto boundary = Centroid::find_boundary(peak_points);
        auto sort_points = [](const Centroid::Point &p1,
                              const Centroid::Point &p2) -> bool {
            return (p1.j < p2.j) || ((p1.j == p2.j) && (p1.i < p2.i));
        };
        std::stable_sort(boundary.begin(), boundary.end(), sort_points);
        CHECK(boundary.size() == expected_boundary[i].size());
        for (size_t j = 0; j < boundary.size(); ++j) {
            CHECK(boundary[j].i == expected_boundary[i][j].i);
            CHECK(boundary[j].j == expected_boundary[i][j].j);
        }
    }
}

TEST_CASE("Find peaks") {
    Grid::Parameters grid_params = {
        {}, {0, 100, 0, 100}, {50.0, 5.0, 5.0}, Instrument::QUAD, 0x00};
    Grid::calculate_dimensions(grid_params);
    std::vector<double> data(grid_params.dimensions.n *
                             grid_params.dimensions.m);
    // Simulating one isotope with three adjacent peaks and an artifact peak.
    // Peak 1
    CHECK(Grid::splat({20, 40, 3}, grid_params, data));
    CHECK(Grid::splat({20, 50, 5}, grid_params, data));
    CHECK(Grid::splat({20, 57, 5.5}, grid_params, data));
    CHECK(Grid::splat({20, 60, 6}, grid_params, data));
    CHECK(Grid::splat({20, 63, 5}, grid_params, data));
    CHECK(Grid::splat({20, 65, 4}, grid_params, data));
    // Peak 2
    CHECK(Grid::splat({50, 50, 4}, grid_params, data));
    CHECK(Grid::splat({50, 57, 4.5}, grid_params, data));
    CHECK(Grid::splat({50, 60, 6}, grid_params, data));
    CHECK(Grid::splat({50, 63, 4}, grid_params, data));
    CHECK(Grid::splat({50, 65, 3}, grid_params, data));
    // Peak 3
    CHECK(Grid::splat({80, 57, 4.5}, grid_params, data));
    CHECK(Grid::splat({80, 60, 5}, grid_params, data));
    CHECK(Grid::splat({80, 63, 3}, grid_params, data));
    CHECK(Grid::splat({80, 65, 2}, grid_params, data));
    // Artifact peak
    CHECK(Grid::splat({10, 10, 4}, grid_params, data));

    auto local_max_points =
        Centroid::find_local_maxima({0, 0, grid_params}, data);
    CHECK(local_max_points.size() == 4);

    std::vector<size_t> expected_size = {
        50,
        40,
        35,
        25,
    };

    std::vector<Centroid::Peak> expected_peaks = {
        {
            4,         // i
            12,        // j
            20,        // mz
            60,        // rt
            17.8731,   // height
            175.126,   // total_intensity
            4.80706,   // sigma_mz
            8.85391,   // sigma_rt
            19.4237,   // mz_centroid
            57.4271,   // rt_centroid
            8.76881,   // height_centroid
            160.357,   // total_intensity_centroid
            0.814667,  // border_background
            {},        // points
            {},        // boundary
        },
        {
            10,        // i
            12,        // j
            50,        // mz
            60,        // rt
            15.4607,   // height
            132.077,   // total_intensity
            4.80706,   // sigma_mz
            6.82431,   // sigma_rt
            49.4237,   // mz_centroid
            58.7364,   // rt_centroid
            8.51286,   // height_centroid
            122.656,   // total_intensity_centroid
            0.749546,  // border_background
            {},        // points
            {},        // boundary
        },
        {
            16,        // i
            12,        // j
            80,        // mz
            60,        // rt
            12.4776,   // height
            88.9893,   // total_intensity
            4.80706,   // sigma_mz
            5.53369,   // sigma_rt
            79.4237,   // mz_centroid
            60.1546,   // rt_centroid
            6.90168,   // height_centroid
            83.5048,   // total_intensity_centroid
            0.542793,  // border_background
            {},        // points
            {},        // boundary
        },
        {
            2,         // i
            2,         // j
            10,        // mz
            10,        // rt
            4,         // height
            24.6757,   // total_intensity
            4.80706,   // sigma_mz
            4.80706,   // sigma_rt
            9.4237,    // mz_centroid
            9.4237,    // rt_centroid
            2.28256,   // height_centroid
            22.0599,   // total_intensity_centroid
            0.317821,  // border_background
            {},        // points
            {},        // boundary
        },
    };

    auto round_double = [](double d) {
        return (long long int)(d * 1000.0) / 1000.0;
    };
    for (size_t i = 0; i < local_max_points.size(); ++i) {
        auto local_max = local_max_points[i];
        auto peak = Centroid::build_peak(local_max, {0, 0, grid_params}, data);
        CHECK(peak.i == expected_peaks[i].i);
        CHECK(peak.j == expected_peaks[i].j);
        CHECK(round_double(peak.mz) == round_double(expected_peaks[i].mz));
        CHECK(round_double(peak.rt) == round_double(expected_peaks[i].rt));
        CHECK(round_double(peak.height) ==
              round_double(expected_peaks[i].height));
        CHECK(round_double(peak.total_intensity) ==
              round_double(expected_peaks[i].total_intensity));
        CHECK(round_double(peak.sigma_mz) ==
              round_double(expected_peaks[i].sigma_mz));
        CHECK(round_double(peak.sigma_rt) ==
              round_double(expected_peaks[i].sigma_rt));
        CHECK(round_double(peak.mz_centroid) ==
              round_double(expected_peaks[i].mz_centroid));
        CHECK(round_double(peak.rt_centroid) ==
              round_double(expected_peaks[i].rt_centroid));
        CHECK(round_double(peak.height_centroid) ==
              round_double(expected_peaks[i].height_centroid));
        CHECK(round_double(peak.total_intensity_centroid) ==
              round_double(expected_peaks[i].total_intensity_centroid));
        CHECK(round_double(peak.border_background) ==
              round_double(expected_peaks[i].border_background));
    }
}
