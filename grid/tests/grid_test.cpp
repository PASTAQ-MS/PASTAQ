#include <string>
#include <vector>

#include "doctest.h"
#include "grid/grid.hpp"

TEST_CASE("Gaussian splatting") {
    // For testing purposes we will set the bounds to [-3 * sigma, +3 * sigma]
    // on both dimensions, even if there is no such thing as negative retention
    // time or mass.
    double sigma_mz = 1.0;
    double sigma_rt = 1.0;

    SUBCASE("Splat on the cencer of the mesh") {
        Grid::Parameters parameters = {
            {},
            {-3.0 * sigma_rt, 3.0 * sigma_rt, -3.0 * sigma_mz, 3.0 * sigma_mz},
            {1.0, sigma_mz, 1.0},
            Instrument::QUAD,
            0x00};
        CHECK(calculate_dimensions(parameters));
        std::vector<double> data(parameters.dimensions.n *
                                 parameters.dimensions.m);
        Grid::splat({0.0, 0.0, 1}, parameters, data);
        std::vector<std::string> expected_weights = {
            "0.000000", "0.000000", "0.000000", "0.000000",
            "0.000000", "0.000000", "0.000000",  // Row 0
            "0.000000", "0.018316", "0.082085", "0.135335",
            "0.082085", "0.018316", "0.000000",  // Row 1
            "0.000000", "0.082085", "0.367879", "0.606531",
            "0.367879", "0.082085", "0.000000",  // Row 2
            "0.000000", "0.135335", "0.606531", "1.000000",
            "0.606531", "0.135335", "0.000000",  // Row 3
            "0.000000", "0.082085", "0.367879", "0.606531",
            "0.367879", "0.082085", "0.000000",  // Row 4
            "0.000000", "0.018316", "0.082085", "0.135335",
            "0.082085", "0.018316", "0.000000",  // Row 5
            "0.000000", "0.000000", "0.000000", "0.000000",
            "0.000000", "0.000000", "0.000000",  // Row 6
        };
        for (size_t j = 0; j < 7; j++) {
            for (size_t i = 0; i < 7; ++i) {
                CHECK(expected_weights[i + 7 * j] ==
                      std::to_string(data[i + 7 * j]));
            }
        }
    }

    SUBCASE("Splatting outside the grid fails") {
        Grid::Parameters parameters = {
            {7, 7},
            {-3.0 * sigma_rt, 3.0 * sigma_rt, -3.0 * sigma_mz, 3.0 * sigma_mz},
            {1.0, sigma_mz, 1.0},
            Instrument::QUAD,
            0x00};
        std::vector<double> data(parameters.dimensions.n *
                                 parameters.dimensions.m);
        CHECK_FALSE(Grid::splat({-20.0, -20.0, 1}, parameters, data));
        CHECK_FALSE(Grid::splat({20.0, 20.0, 1}, parameters, data));
        CHECK_FALSE(Grid::splat({-20.0, 20.0, 1}, parameters, data));
        CHECK_FALSE(Grid::splat({20.0, +20.0, 1}, parameters, data));
    }

    SUBCASE("Splat outside but range falls inside") {
        Grid::Parameters parameters = {
            {7, 7},
            {-3.0 * sigma_rt, 3.0 * sigma_rt, -3.0 * sigma_mz, 3.0 * sigma_mz},
            {1.0, sigma_mz, 1.0},
            Instrument::QUAD,
            0x00};
        std::vector<double> data(parameters.dimensions.n *
                                 parameters.dimensions.m);
        CHECK(Grid::splat({-5.0, 0.0, 1}, parameters, data));
        CHECK(Grid::splat({5.0, 0.0, 1}, parameters, data));
        CHECK(Grid::splat({0.0, -5.0, 1}, parameters, data));
        CHECK(Grid::splat({0.0, 5.0, 1}, parameters, data));
        std::vector<std::string> expected_weights = {
            "0.000000", "0.018316", "0.082085", "0.135335",
            "0.082085", "0.018316", "0.000000",  // Row 0
            "0.018316", "0.000000", "0.000000", "0.000000",
            "0.000000", "0.000000", "0.018316",  // Row 1
            "0.082085", "0.000000", "0.000000", "0.000000",
            "0.000000", "0.000000", "0.082085",  // Row 2
            "0.135335", "0.000000", "0.000000", "0.000000",
            "0.000000", "0.000000", "0.135335",  // Row 3
            "0.082085", "0.000000", "0.000000", "0.000000",
            "0.000000", "0.000000", "0.082085",  // Row 4
            "0.018316", "0.000000", "0.000000", "0.000000",
            "0.000000", "0.000000", "0.018316",  // Row 5
            "0.000000", "0.018316", "0.082085", "0.135335",
            "0.082085", "0.018316", "0.000000",  // Row 6
        };
        for (size_t j = 0; j < 7; j++) {
            for (size_t i = 0; i < 7; ++i) {
                CHECK(expected_weights[i + 7 * j] ==
                      std::to_string(data[i + 7 * j]));
            }
        }
    }
}

TEST_CASE("Test the fetching of real world parameters from index") {
    SUBCASE("Check for proper interpolation") {
        Grid::Parameters parameters = {{4, 4},
                                       {0.0, 75.0, 200.0, 800.0},
                                       {200, 0.01, 1.0},
                                       Instrument::QUAD,
                                       0x00};
        CHECK(Grid::mz_at(0, parameters) == 200.0);
        CHECK(Grid::mz_at(1, parameters) == 400.0);
        CHECK(Grid::mz_at(2, parameters) == 600.0);
        CHECK(Grid::mz_at(3, parameters) == 800.0);
        CHECK(Grid::rt_at(0, parameters) == 0.0);
        CHECK(Grid::rt_at(1, parameters) == 25.0);
        CHECK(Grid::rt_at(2, parameters) == 50.0);
        CHECK(Grid::rt_at(3, parameters) == 75.0);
    };
    SUBCASE("Using a warped grid") {
        auto round_double = [](double d) {
            return (long long int)(d * 1000.0) / 1000.0;
        };
        SUBCASE("QUAD") {
            Grid::Parameters parameters = {{4, 4},
                                           {0.0, 75.0, 200.0, 800.0},
                                           {200, 1.0, 1.0},
                                           Instrument::QUAD,
                                           Grid::Flags::WARPED_MESH};
            CHECK(Grid::mz_at(0, parameters) == 200.0);
            CHECK(Grid::mz_at(1, parameters) == 400.0);
            CHECK(Grid::mz_at(2, parameters) == 600.0);
            CHECK(Grid::mz_at(3, parameters) == 800.0);
        }
        SUBCASE("ORBITRAP") {
            Grid::Parameters parameters = {{},
                                           {0.0, 75.0, 200.0, 800.0},
                                           {200, 1.0, 1.0},
                                           Instrument::ORBITRAP,
                                           Grid::Flags::WARPED_MESH};
            Grid::calculate_dimensions(parameters);
            CHECK(round_double(Grid::mz_at(0, parameters)) == 200.0);
            CHECK(round_double(Grid::mz_at(1, parameters)) == 202.015);
            CHECK(round_double(Grid::mz_at(2, parameters)) == 204.06);
            CHECK(round_double(Grid::mz_at(3, parameters)) == 206.137);
            CHECK(round_double(Grid::mz_at(100, parameters)) == 800.0);
        }
        SUBCASE("FTICR") {
            Grid::Parameters parameters = {{},
                                           {0.0, 75.0, 200.0, 800.0},
                                           {200, 1.0, 1.0},
                                           Instrument::FTICR,
                                           Grid::Flags::WARPED_MESH};
            Grid::calculate_dimensions(parameters);
            CHECK(round_double(Grid::mz_at(0, parameters)) == 200.0);
            CHECK(round_double(Grid::mz_at(1, parameters)) == 201.005);
            CHECK(round_double(Grid::mz_at(2, parameters)) == 202.02);
            CHECK(round_double(Grid::mz_at(3, parameters)) == 203.045);
            CHECK(round_double(Grid::mz_at(100, parameters)) == 400);
            CHECK(round_double(Grid::mz_at(150, parameters)) == 800);
        }
        SUBCASE("TOF") {
            Grid::Parameters parameters = {{},
                                           {0.0, 75.0, 200.0, 800.0},
                                           {200, 1.0, 1.0},
                                           Instrument::TOF,
                                           Grid::Flags::WARPED_MESH};
            Grid::calculate_dimensions(parameters);
            CHECK(round_double(Grid::mz_at(0, parameters)) == 200.0);
            CHECK(round_double(Grid::mz_at(1, parameters)) == 201.002);
            CHECK(round_double(Grid::mz_at(2, parameters)) == 202.01);
            CHECK(round_double(Grid::mz_at(3, parameters)) == 203.022);
        }
    }
}
