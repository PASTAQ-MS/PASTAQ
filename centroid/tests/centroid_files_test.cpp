#include <iostream>

#include "doctest.h"

#include "centroid/centroid_files.hpp"
#include "testutils/mock_stream.hpp"

TEST_CASE("Read/Write header from/to stream") {
    std::vector<Centroid::Files::Bpks::Header> source_data = {
        {1, 2, 3, {}},
        {4, 5, 6, {}},
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
        CHECK(read_value.header_length == value.header_length);
        CHECK(read_value.num_peaks == value.num_peaks);
        CHECK(read_value.spec_version == value.spec_version);
        // TODO(alex): Verify Header::Grid::Parameters
    }
}
