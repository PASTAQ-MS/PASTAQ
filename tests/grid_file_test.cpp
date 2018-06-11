#include <fstream>
#include <iostream>
#include <streambuf>

#include "doctest.h"
#include "grid_file.hpp"
#include "utils/mock_stream.hpp"

TEST_CASE("Writing functions to stream") {
    SUBCASE("uint32") {
        std::vector<unsigned int> data_source = {1, 2, 3, 4, 5, 6, 7, 8, 9, 0};
        std::vector<unsigned int> data_destination(data_source.size());
        MockStream<unsigned int> stream(data_destination);
        for (const auto& e : data_source) {
            Grid::File::save_uint32(stream, e);
        }
        for (int i = 0; i < data_source.size(); ++i) {
            unsigned int j = 0;
            Grid::File::load_uint32(stream, &j);
            CHECK(data_destination[i] == data_source[i]);
            CHECK(j == data_source[i]);
        }
    }
    SUBCASE("double") {
        std::vector<double> data_source = {1, 2, 3, 4, 5, 6, 7, 8, 9, 0};
        std::vector<double> data_destination(data_source.size());
        MockStream<double> stream(data_destination);
        for (const auto& e : data_source) {
            Grid::File::save_double(stream, e);
        }
        for (int i = 0; i < data_source.size(); ++i) {
            double j = 0;
            Grid::File::load_double(stream, &j);
            CHECK(data_destination[i] == data_source[i]);
            CHECK(j == data_source[i]);
        }
    }
}

TEST_CASE("Writing parameters to the stream") {
    SUBCASE("Testing on MockStream") {
        std::vector<char> data(sizeof(Grid::Parameters) +
                               sizeof(Grid::File::Parameters));
        MockStream<char> stream(data);
        Grid::Parameters source_parameters = {
            {6, 2}, {3, 4, 5, 6}, {7, 8, 9}, Instrument::QUAD};
        CHECK(Grid::File::write_parameters(stream, source_parameters));
        Grid::Parameters dest_parameters = {};
        CHECK(Grid::File::load_parameters(stream, &dest_parameters));
        CHECK(dest_parameters.bounds.min_rt == source_parameters.bounds.min_rt);
        CHECK(dest_parameters.bounds.max_rt == source_parameters.bounds.max_rt);
        CHECK(dest_parameters.bounds.min_mz == source_parameters.bounds.min_mz);
        CHECK(dest_parameters.bounds.max_mz == source_parameters.bounds.max_mz);
        CHECK(dest_parameters.dimensions.n == source_parameters.dimensions.n);
        CHECK(dest_parameters.dimensions.m == source_parameters.dimensions.m);
        CHECK(dest_parameters.smoothing_params.mz ==
              source_parameters.smoothing_params.mz);
        CHECK(dest_parameters.smoothing_params.sigma_mz ==
              source_parameters.smoothing_params.sigma_mz);
        CHECK(dest_parameters.smoothing_params.sigma_rt ==
              source_parameters.smoothing_params.sigma_rt);
        CHECK(dest_parameters.instrument_type ==
              source_parameters.instrument_type);
        CHECK(dest_parameters.flags == source_parameters.flags);
    }

    SUBCASE("Testing on file stream") {
        Grid::Parameters source_parameters = {
            {6, 2}, {3, 4, 5, 6}, {7, 8, 9}, Instrument::QUAD};
        std::ofstream fileout("test_dat_file.dat",
                              std::ios::out | std::ios::binary);
        CHECK(Grid::File::write_parameters(fileout, source_parameters));
        fileout.close();
        std::ifstream filein("test_dat_file.dat",
                             std::ios::out | std::ios::binary);
        Grid::Parameters dest_parameters = {};
        CHECK(Grid::File::load_parameters(filein, &dest_parameters));
        CHECK(dest_parameters.bounds.min_rt == source_parameters.bounds.min_rt);
        CHECK(dest_parameters.bounds.max_rt == source_parameters.bounds.max_rt);
        CHECK(dest_parameters.bounds.min_mz == source_parameters.bounds.min_mz);
        CHECK(dest_parameters.bounds.max_mz == source_parameters.bounds.max_mz);
        CHECK(dest_parameters.dimensions.n == source_parameters.dimensions.n);
        CHECK(dest_parameters.dimensions.m == source_parameters.dimensions.m);
        CHECK(dest_parameters.smoothing_params.mz ==
              source_parameters.smoothing_params.mz);
        CHECK(dest_parameters.smoothing_params.sigma_mz ==
              source_parameters.smoothing_params.sigma_mz);
        CHECK(dest_parameters.smoothing_params.sigma_rt ==
              source_parameters.smoothing_params.sigma_rt);
        CHECK(dest_parameters.instrument_type ==
              source_parameters.instrument_type);
        CHECK(dest_parameters.flags == source_parameters.flags);
    }
}

TEST_CASE("Writing data/parameters to the stream") {
    SUBCASE("Testing on MockStream") {
        std::vector<double> source_data = {
            1, 2, 3, 4, 5,  // Row 1
            6, 5, 4, 3, 2,  // Row 2
        };
        Grid::Parameters source_parameters = {
            {5, 2}, {3, 4, 5, 6}, {7, 8, 9}, Instrument::QUAD};
        std::vector<char> stream_data(sizeof(double) * source_data.size() +
                                      sizeof(Grid::Parameters) +
                                      sizeof(Grid::File::Parameters));
        MockStream<char> stream(stream_data);
        CHECK(Grid::File::write(stream, source_data, source_parameters));
        std::vector<double> dest_data = {};
        Grid::Parameters dest_parameters = {};
        CHECK(Grid::File::load(stream, &dest_data, &dest_parameters));

        // Check that the data is restored succesfully.
        for (int i = 0; i < dest_data.size(); ++i) {
            CHECK(source_data[i] == dest_data[i]);
        }

        // Check that the parameters are loaded correctly from the stream.
        CHECK(dest_parameters.bounds.min_rt == source_parameters.bounds.min_rt);
        CHECK(dest_parameters.bounds.max_rt == source_parameters.bounds.max_rt);
        CHECK(dest_parameters.bounds.min_mz == source_parameters.bounds.min_mz);
        CHECK(dest_parameters.bounds.max_mz == source_parameters.bounds.max_mz);
        CHECK(dest_parameters.dimensions.n == source_parameters.dimensions.n);
        CHECK(dest_parameters.dimensions.m == source_parameters.dimensions.m);
        CHECK(dest_parameters.smoothing_params.mz ==
              source_parameters.smoothing_params.mz);
        CHECK(dest_parameters.smoothing_params.sigma_mz ==
              source_parameters.smoothing_params.sigma_mz);
        CHECK(dest_parameters.smoothing_params.sigma_rt ==
              source_parameters.smoothing_params.sigma_rt);
        CHECK(dest_parameters.instrument_type ==
              source_parameters.instrument_type);
        CHECK(dest_parameters.flags == source_parameters.flags);
    }

    SUBCASE("Testing on file stream") {
        std::vector<double> source_data = {
            1, 2, 3, 4, 5,  // Row 1
            6, 5, 4, 3, 2,  // Row 2
        };
        Grid::Parameters source_parameters = {
            {5, 2}, {3, 4, 5, 6}, {7, 8, 9}, Instrument::QUAD};
        std::ofstream fileout("test_dat_file_all.dat",
                              std::ios::out | std::ios::binary);
        CHECK(Grid::File::write(fileout, source_data, source_parameters));
        fileout.close();
        std::ifstream filein("test_dat_file_all.dat",
                             std::ios::in | std::ios::binary);
        std::vector<double> dest_data = {};
        Grid::Parameters dest_parameters = {};
        CHECK(Grid::File::load(filein, &dest_data, &dest_parameters));

        // Check that the data is restored succesfully.
        for (int i = 0; i < dest_data.size(); ++i) {
            CHECK(source_data[i] == dest_data[i]);
        }

        // Check that the parameters are loaded correctly from the stream.
        CHECK(dest_parameters.bounds.min_rt == source_parameters.bounds.min_rt);
        CHECK(dest_parameters.bounds.max_rt == source_parameters.bounds.max_rt);
        CHECK(dest_parameters.bounds.min_mz == source_parameters.bounds.min_mz);
        CHECK(dest_parameters.bounds.max_mz == source_parameters.bounds.max_mz);
        CHECK(dest_parameters.dimensions.n == source_parameters.dimensions.n);
        CHECK(dest_parameters.dimensions.m == source_parameters.dimensions.m);
        CHECK(dest_parameters.smoothing_params.mz ==
              source_parameters.smoothing_params.mz);
        CHECK(dest_parameters.smoothing_params.sigma_mz ==
              source_parameters.smoothing_params.sigma_mz);
        CHECK(dest_parameters.smoothing_params.sigma_rt ==
              source_parameters.smoothing_params.sigma_rt);
        CHECK(dest_parameters.instrument_type ==
              source_parameters.instrument_type);
        CHECK(dest_parameters.flags == source_parameters.flags);
    }
}
