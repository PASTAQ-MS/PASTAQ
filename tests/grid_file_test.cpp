#include <fstream>
#include <iostream>
#include <streambuf>

#include "doctest.h"
#include "grid_file.hpp"
#include "utils/mock_stream.hpp"

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

    SUBCASE("Loading custom range") {
        std::vector<double> source_data = {
            1, 2, 3, 4, 5,  // Row 1
            6, 5, 4, 3, 2,  // Row 2
        };
        Grid::Parameters source_parameters = {
            {5, 2}, {0, 1, 0, 4}, {1, 1, 1}, Instrument::QUAD};
        std::vector<char> stream_data(sizeof(double) * source_data.size() +
                                      sizeof(Grid::Parameters) +
                                      sizeof(Grid::File::Parameters));
        MockStream<char> stream(stream_data);
        CHECK(Grid::File::write(stream, source_data, source_parameters));
        std::vector<double> dest_data = {};
        Grid::Parameters dest_parameters = {};
        auto load_range_results = Grid::File::load_range(
            stream, {0, 1, 1, 3}, &dest_data, &dest_parameters);
        CHECK(load_range_results);
        if (load_range_results) {
            std::vector<double> expected = {
                2, 3, 4,  // Row 1
                5, 4, 3,  // Row 2
            };
            for (size_t i = 0; i < dest_data.size(); ++i) {
                CHECK(expected[i] == dest_data[i]);
            }
            CHECK(dest_parameters.dimensions.n == 3);
            CHECK(dest_parameters.dimensions.m == 2);
            CHECK(dest_parameters.bounds.min_mz == 1);
            CHECK(dest_parameters.bounds.max_mz == 3);
            CHECK(dest_parameters.bounds.min_rt == 0);
            CHECK(dest_parameters.bounds.max_rt == 1);
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
