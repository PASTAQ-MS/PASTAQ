#include <iostream>

#include "centroid/centroid_files.hpp"
#include "doctest.h"
#include "testutils/mock_stream.hpp"

TEST_CASE("DUMMY") {
    // Grid::Parameters parameters = {
    //{}, {0, 100, 0, 100}, {50.0, 5.0, 5.0}, Instrument::QUAD, 0x00};
    // Grid::calculate_dimensions(parameters);
    // std::vector<double> data(parameters.dimensions.n *
    // parameters.dimensions.m);
    //// Simulating one isotope with three adjacent peaks and an artifact peak.
    //// Peak 1
    // CHECK(Grid::splat({20, 40, 3}, parameters, data));
    // CHECK(Grid::splat({20, 50, 5}, parameters, data));
    // CHECK(Grid::splat({20, 57, 5.5}, parameters, data));
    // CHECK(Grid::splat({20, 60, 6}, parameters, data));
    // CHECK(Grid::splat({20, 63, 5}, parameters, data));
    // CHECK(Grid::splat({20, 65, 4}, parameters, data));
    //// Peak 2
    // CHECK(Grid::splat({50, 50, 4}, parameters, data));
    // CHECK(Grid::splat({50, 57, 4.5}, parameters, data));
    // CHECK(Grid::splat({50, 60, 6}, parameters, data));
    // CHECK(Grid::splat({50, 63, 4}, parameters, data));
    // CHECK(Grid::splat({50, 65, 3}, parameters, data));
    //// Peak 3
    // CHECK(Grid::splat({80, 57, 4.5}, parameters, data));
    // CHECK(Grid::splat({80, 60, 5}, parameters, data));
    // CHECK(Grid::splat({80, 63, 3}, parameters, data));
    // CHECK(Grid::splat({80, 65, 2}, parameters, data));
    //// Artifact peak
    // CHECK(Grid::splat({10, 10, 4}, parameters, data));

    // auto points = Centroid::find_local_maxima(parameters, data);
    // CHECK(points.size() == 4);
    // CHECK(1 == 0);
}

// TODO(alex): move to file utils test module.
TEST_CASE("Read/Write double values") {
    std::vector<double> source_data = {1.0, 2.0, 3.0, 4.0, -1.0, -2.0, -3.0};
    std::vector<double> destination_data;
    std::vector<char> stream_data(sizeof(double) * source_data.size());
    MockStream<char> stream(stream_data);
    // Write data to stream.
    for (const auto &value : source_data) {
        CHECK(Centroid::Files::Bpks::write_double(stream, value));
    }
    // Read data from stream.
    for (const auto &value : source_data) {
        double read_value;
        CHECK(Centroid::Files::Bpks::read_double(stream, &read_value));
        CHECK(read_value == value);
    }
}
