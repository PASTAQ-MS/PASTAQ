#include <iomanip>
#include <iostream>

#include "centroid.hpp"
#include "doctest.h"

TEST_CASE("Find local maxima") {
    Grid::Parameters parameters = {
        {}, {0, 100, 0, 100}, {50.0, 5.0, 5.0}, Instrument::QUAD, 0x00};
    Grid::calculate_dimensions(parameters);
    std::vector<double> data(parameters.dimensions.n * parameters.dimensions.m);
    // Simulating one isotope with three adjacent peaks and an artifact peak.
    // Peak 1
    CHECK(Grid::splat({20, 40, 3}, parameters, data));
    CHECK(Grid::splat({20, 50, 5}, parameters, data));
    CHECK(Grid::splat({20, 57, 5.5}, parameters, data));
    CHECK(Grid::splat({20, 60, 6}, parameters, data));
    CHECK(Grid::splat({20, 63, 5}, parameters, data));
    CHECK(Grid::splat({20, 65, 4}, parameters, data));
    // Peak 2
    CHECK(Grid::splat({50, 50, 4}, parameters, data));
    CHECK(Grid::splat({50, 57, 4.5}, parameters, data));
    CHECK(Grid::splat({50, 60, 6}, parameters, data));
    CHECK(Grid::splat({50, 63, 4}, parameters, data));
    CHECK(Grid::splat({50, 65, 3}, parameters, data));
    // Peak 3
    CHECK(Grid::splat({80, 57, 4.5}, parameters, data));
    CHECK(Grid::splat({80, 60, 5}, parameters, data));
    CHECK(Grid::splat({80, 63, 3}, parameters, data));
    CHECK(Grid::splat({80, 65, 2}, parameters, data));
    // Artifact peak
    CHECK(Grid::splat({10, 10, 4}, parameters, data));

    auto points = Centroid::find_local_maxima(parameters, data);
    CHECK(points.size() == 4);
}

TEST_CASE("Find peak points") {
    Grid::Parameters parameters = {
        {}, {0, 100, 0, 100}, {50.0, 5.0, 5.0}, Instrument::QUAD, 0x00};
    Grid::calculate_dimensions(parameters);
    std::vector<double> data(parameters.dimensions.n * parameters.dimensions.m);
    // Simulating one isotope with three adjacent peaks and an artifact peak.
    // Peak 1
    CHECK(Grid::splat({20, 40, 3}, parameters, data));
    CHECK(Grid::splat({20, 50, 5}, parameters, data));
    CHECK(Grid::splat({20, 57, 5.5}, parameters, data));
    CHECK(Grid::splat({20, 60, 6}, parameters, data));
    CHECK(Grid::splat({20, 63, 5}, parameters, data));
    CHECK(Grid::splat({20, 65, 4}, parameters, data));
    // Peak 2
    CHECK(Grid::splat({50, 50, 4}, parameters, data));
    CHECK(Grid::splat({50, 57, 4.5}, parameters, data));
    CHECK(Grid::splat({50, 60, 6}, parameters, data));
    CHECK(Grid::splat({50, 63, 4}, parameters, data));
    CHECK(Grid::splat({50, 65, 3}, parameters, data));
    // Peak 3
    CHECK(Grid::splat({80, 57, 4.5}, parameters, data));
    CHECK(Grid::splat({80, 60, 5}, parameters, data));
    CHECK(Grid::splat({80, 63, 3}, parameters, data));
    CHECK(Grid::splat({80, 65, 2}, parameters, data));
    // Artifact peak
    CHECK(Grid::splat({10, 10, 4}, parameters, data));

    auto local_max_points = Centroid::find_local_maxima(parameters, data);
    CHECK(local_max_points.size() == 4);

    std::vector<size_t> expected_size = {
        25,
        50,
        40,
        35,
    };

    for (size_t i = 0; i < local_max_points.size(); ++i) {
        auto local_max = local_max_points[i];
        auto peak_points =
            Centroid::find_peak_points(local_max, parameters, data);
        CHECK(expected_size[i] == peak_points.size());
    }
}

TEST_CASE("Find peak boundaries") {}
