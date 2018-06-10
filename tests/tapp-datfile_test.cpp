#include <iostream>
#include <streambuf>

#include "doctest.h"
#include "tapp-datfile.hpp"
#include "utils/mock_stream.hpp"

TEST_CASE("Writing functions to stream") {
    SUBCASE("uint32") {
        std::vector<unsigned int> data_source = {1, 2, 3, 4, 5, 6, 7, 8, 9, 0};
        std::vector<unsigned int> data_destination(data_source.size());
        MockStream<unsigned int> stream(data_destination);
        for (const auto& e : data_source) {
            DatFile::save_uint32(stream, e);
        }
        for (int i = 0; i < data_source.size(); ++i) {
            unsigned int j = 0;
            DatFile::load_uint32(stream, &j);
            CHECK(data_destination[i] == data_source[i]);
            CHECK(j == data_source[i]);
        }
    }
    SUBCASE("double") {
        std::vector<double> data_source = {1, 2, 3, 4, 5, 6, 7, 8, 9, 0};
        std::vector<double> data_destination(data_source.size());
        MockStream<double> stream(data_destination);
        for (const auto& e : data_source) {
            DatFile::save_double(stream, e);
        }
        for (int i = 0; i < data_source.size(); ++i) {
            double j = 0;
            DatFile::load_double(stream, &j);
            CHECK(data_destination[i] == data_source[i]);
            CHECK(j == data_source[i]);
        }
    }
}

TEST_CASE("Writing parameters to the stream") {
    std::vector<char> data(255);
    MockStream<char> stream(data);
    Grid::Parameters parameters = {
        {1, 2}, {3, 4, 5, 6}, {7, 8, 9}, Instrument::QUAD};
    CHECK(DatFile::write_parameters(stream, parameters));
    //CHECK(0 == 1);
}
