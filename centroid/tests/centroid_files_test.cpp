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

TEST_CASE("Read/Write found peaks to bpks stream") {
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

    std::vector<Centroid::Peak> expected_peaks = {
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
    };

    // TODO(alex): These should be generalized to testutils module
    auto round_double = [](double d) {
        return (long long int)(d * 1000.0) / 1000.0;
    };
    std::vector<Centroid::Peak> peaks;
    for (size_t i = 0; i < local_max_points.size(); ++i) {
        auto local_max = local_max_points[i];
        auto peak = Centroid::build_peak(local_max, parameters, data);
        peaks.push_back(peak);
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
    // Writing data to stream.
    std::vector<char> stream_data(1024 * 1024);
    MockStream<char> stream(stream_data);
    CHECK(Centroid::Files::Bpks::write_peaks(stream, parameters, peaks));
    // Reading data from stream.
    std::vector<Centroid::Peak> read_peaks;
    Grid::Parameters read_parameters;
    CHECK(Centroid::Files::Bpks::read_peaks(stream, &read_parameters,
                                            &read_peaks));

    // Check that parameters were read successfully.
    CHECK(parameters.dimensions.n == read_parameters.dimensions.n);
    CHECK(parameters.dimensions.m == read_parameters.dimensions.m);
    CHECK(round_double(parameters.bounds.min_rt) ==
          round_double(read_parameters.bounds.min_rt));
    CHECK(round_double(parameters.bounds.max_rt) ==
          round_double(read_parameters.bounds.max_rt));
    CHECK(round_double(parameters.bounds.min_mz) ==
          round_double(read_parameters.bounds.min_mz));
    CHECK(round_double(parameters.bounds.max_mz) ==
          round_double(read_parameters.bounds.max_mz));
    CHECK(round_double(parameters.smoothing_params.mz) ==
          round_double(read_parameters.smoothing_params.mz));
    CHECK(round_double(parameters.smoothing_params.sigma_mz) ==
          round_double(read_parameters.smoothing_params.sigma_mz));
    CHECK(round_double(parameters.smoothing_params.sigma_rt) ==
          round_double(read_parameters.smoothing_params.sigma_rt));
    CHECK(parameters.instrument_type == read_parameters.instrument_type);
    CHECK(parameters.flags == read_parameters.flags);

    // Check that peaks were read successfully.
    CHECK(peaks.size() == read_peaks.size());
    for (size_t i = 0; i < peaks.size(); ++i) {
        CHECK(peaks[i].i == read_peaks[i].i);
        CHECK(peaks[i].j == read_peaks[i].j);
        CHECK(round_double(peaks[i].mz) == round_double(read_peaks[i].mz));
        CHECK(round_double(peaks[i].rt) == round_double(read_peaks[i].rt));
        CHECK(round_double(peaks[i].total_intensity) ==
              round_double(read_peaks[i].total_intensity));
        CHECK(round_double(peaks[i].sigma_mz) ==
              round_double(read_peaks[i].sigma_mz));
        CHECK(round_double(peaks[i].sigma_rt) ==
              round_double(read_peaks[i].sigma_rt));
        CHECK(round_double(peaks[i].mz_centroid) ==
              round_double(read_peaks[i].mz_centroid));
        CHECK(round_double(peaks[i].rt_centroid) ==
              round_double(read_peaks[i].rt_centroid));
        CHECK(round_double(peaks[i].height_centroid) ==
              round_double(read_peaks[i].height_centroid));
        CHECK(round_double(peaks[i].total_intensity_centroid) ==
              round_double(read_peaks[i].total_intensity_centroid));
        CHECK(round_double(peaks[i].border_background) ==
              round_double(read_peaks[i].border_background));

        CHECK(peaks[i].points.size() == read_peaks[i].points.size());
        for (size_t j = 0; j < peaks[i].points.size(); ++j) {
            CHECK(peaks[i].points[j].i == read_peaks[i].points[j].i);
            CHECK(peaks[i].points[j].j == read_peaks[i].points[j].j);
            CHECK(round_double(peaks[i].points[j].height) ==
                  round_double(read_peaks[i].points[j].height));
        }

        CHECK(peaks[i].boundary.size() == read_peaks[i].boundary.size());
        for (size_t j = 0; j < peaks[i].boundary.size(); ++j) {
            CHECK(peaks[i].boundary[j].i == read_peaks[i].boundary[j].i);
            CHECK(peaks[i].boundary[j].j == read_peaks[i].boundary[j].j);
            CHECK(round_double(peaks[i].boundary[j].height) ==
                  round_double(read_peaks[i].boundary[j].height));
        }
    }
}
