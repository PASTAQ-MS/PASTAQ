#include "grid_runners.hpp"
#include "doctest.h"

TEST_CASE("Parallel and serial execution offer the same results") {
    auto round_double = [](double d) {
        return (long long int)(d * 1000.0) / 1000.0;
    };

    Grid::Parameters parameters = {
        {0, 0}, {0.0, 100, 0.0, 100}, {200, 1.0, 1.0}, Instrument::QUAD};
    calculate_dimensions(parameters);

    // Generate dummy peaks for testing.
    std::vector<Grid::Peak> dummy_peaks;
    for (size_t i = 0; i < 100; ++i) {
        double d = i;
        dummy_peaks.push_back({d, d, d});
    }
    auto results_serial = Grid::Runners::Serial::run(parameters, dummy_peaks);
    auto results_parallel =
        Grid::Runners::Parallel::run(3, parameters, dummy_peaks);
    CHECK(results_serial.size() == results_parallel.size());
    for (size_t i = 0; i < results_serial.size(); ++i) {
        INFO(i);
        CHECK(round_double(results_serial[i]) ==
              round_double(results_parallel[i]));
    }
}
