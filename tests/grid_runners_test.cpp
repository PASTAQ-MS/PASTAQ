#include "grid_runners.hpp"
#include "doctest.h"

TEST_CASE("Parallel and serial execution offer the same results") {
    SUBCASE("NON WARPED") {
        auto round_double = [](double d) {
            return (long long int)(d * 1000.0) / 1000.0;
        };

        Grid::Parameters parameters = {{0, 0},
                                       {0.0, 1000, 0.0, 1000},
                                       {200, 1.0, 1.0},
                                       Instrument::ORBITRAP};
        calculate_dimensions(parameters);

        // Generate dummy peaks for testing.
        std::vector<Grid::Point> dummy_peaks;
        for (size_t i = 0; i < 1000; ++i) {
            double d = i;
            dummy_peaks.push_back({d, d, d});
        }
        auto results_serial =
            Grid::Runners::Serial::run(parameters, dummy_peaks);
        // Testing the results stay stable regardless of the number of threads
        // used.
        {
            auto results_parallel =
                Grid::Runners::Parallel::run(1, parameters, dummy_peaks);
            CHECK(results_serial.size() == results_parallel.size());
            for (size_t i = 0; i < results_serial.size(); ++i) {
                CHECK(round_double(results_serial[i]) ==
                      round_double(results_parallel[i]));
            }
        }
        {
            auto results_parallel =
                Grid::Runners::Parallel::run(4, parameters, dummy_peaks);
            CHECK(results_serial.size() == results_parallel.size());
            for (size_t i = 0; i < results_serial.size(); ++i) {
                CHECK(round_double(results_serial[i]) ==
                      round_double(results_parallel[i]));
            }
        }
        {
            auto results_parallel =
                Grid::Runners::Parallel::run(8, parameters, dummy_peaks);
            CHECK(results_serial.size() == results_parallel.size());
            for (size_t i = 0; i < results_serial.size(); ++i) {
                CHECK(round_double(results_serial[i]) ==
                      round_double(results_parallel[i]));
            }
        }
        {
            auto results_parallel =
                Grid::Runners::Parallel::run(12, parameters, dummy_peaks);
            CHECK(results_serial.size() == results_parallel.size());
            for (size_t i = 0; i < results_serial.size(); ++i) {
                CHECK(round_double(results_serial[i]) ==
                      round_double(results_parallel[i]));
            }
        }
    }

    SUBCASE("WARPED") {
        auto round_double = [](double d) {
            return (long long int)(d * 1000.0) / 1000.0;
        };

        Grid::Parameters parameters = {{0, 0},
                                       {0.0, 1000, 0.0, 1000},
                                       {200, 1.0, 1.0},
                                       Instrument::ORBITRAP,
                                       Grid::Flags::WARPED_MESH};
        calculate_dimensions(parameters);

        // Generate dummy peaks for testing.
        std::vector<Grid::Point> dummy_peaks;
        for (size_t i = 0; i < 1000; ++i) {
            double d = i;
            dummy_peaks.push_back({d, d, d});
        }
        auto results_serial =
            Grid::Runners::Serial::run(parameters, dummy_peaks);
        // Testing the results stay stable regardless of the number of threads
        // used.
        {
            auto results_parallel =
                Grid::Runners::Parallel::run(1, parameters, dummy_peaks);
            CHECK(results_serial.size() == results_parallel.size());
            for (size_t i = 0; i < results_serial.size(); ++i) {
                CHECK(round_double(results_serial[i]) ==
                      round_double(results_parallel[i]));
            }
        }
        {
            auto results_parallel =
                Grid::Runners::Parallel::run(4, parameters, dummy_peaks);
            CHECK(results_serial.size() == results_parallel.size());
            for (size_t i = 0; i < results_serial.size(); ++i) {
                CHECK(round_double(results_serial[i]) ==
                      round_double(results_parallel[i]));
            }
        }
        {
            auto results_parallel =
                Grid::Runners::Parallel::run(8, parameters, dummy_peaks);
            CHECK(results_serial.size() == results_parallel.size());
            for (size_t i = 0; i < results_serial.size(); ++i) {
                CHECK(round_double(results_serial[i]) ==
                      round_double(results_parallel[i]));
            }
        }
        {
            auto results_parallel =
                Grid::Runners::Parallel::run(12, parameters, dummy_peaks);
            CHECK(results_serial.size() == results_parallel.size());
            for (size_t i = 0; i < results_serial.size(); ++i) {
                CHECK(round_double(results_serial[i]) ==
                      round_double(results_parallel[i]));
            }
        }
    }
}
