#include <iostream>

#include "doctest.h"

#include "centroid/centroid_files.hpp"
#include "testutils/mock_stream.hpp"

TEST_CASE("Read/Write header from/to stream") {
    std::vector<Centroid::Files::Bpks::Header> source_data = {
        {1, 2, {{1, 2}, {3, 4, 5, 6}, {7, 8, 9}, 0, 1}},
        {4, 5, {{3, 1}, {2, 9, 1, 2}, {8, 0, 7}, 2, 3}},
    };
    std::vector<Centroid::Files::Bpks::Header> destination_data;
    std::vector<char> stream_data(sizeof(Centroid::Files::Bpks::Header) *
                                  source_data.size());
    MockStream<char> stream(stream_data);
    // Write data to stream.
    for (const auto &value : source_data) {
        CHECK(Centroid::Files::Bpks::write_header(stream, value));
    }
    // Read data from stream.
    for (const auto &value : source_data) {
        Centroid::Files::Bpks::Header read_value;
        CHECK(Centroid::Files::Bpks::read_header(stream, &read_value));
        CHECK(read_value.num_peaks == value.num_peaks);
        CHECK(read_value.spec_version == value.spec_version);
        CHECK(read_value.grid_params.dimensions.n ==
              value.grid_params.dimensions.n);
        CHECK(read_value.grid_params.dimensions.m ==
              value.grid_params.dimensions.m);
        CHECK(read_value.grid_params.bounds.min_mz ==
              value.grid_params.bounds.min_mz);
        CHECK(read_value.grid_params.bounds.max_mz ==
              value.grid_params.bounds.max_mz);
        CHECK(read_value.grid_params.bounds.min_rt ==
              value.grid_params.bounds.min_rt);
        CHECK(read_value.grid_params.bounds.max_rt ==
              value.grid_params.bounds.max_rt);
        CHECK(read_value.grid_params.smoothing_params.mz ==
              value.grid_params.smoothing_params.mz);
        CHECK(read_value.grid_params.smoothing_params.sigma_mz ==
              value.grid_params.smoothing_params.sigma_mz);
        CHECK(read_value.grid_params.smoothing_params.sigma_rt ==
              value.grid_params.smoothing_params.sigma_rt);
    }
}
