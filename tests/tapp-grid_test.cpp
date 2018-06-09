#include <string>
#include <vector>

#include "doctest.h"
#include "tapp-grid.hpp"
#include "tapp-mesh.hpp"

TEST_CASE("Gaussian splatting") {
    // For testing purposes we will set the bounds to [-3 * sigma, +3 * sigma]
    // on both dimensions, even if there is no such thing as negative retention
    // time or mass.
    double sigma_mz = 1.0;
    double sigma_rt = 1.0;

    SUBCASE("Splat on the cencer of the mesh") {
        auto mesh = Mesh(
            {7, 7},
            {-3.0 * sigma_rt, 3.0 * sigma_rt, -3.0 * sigma_mz, 3.0 * sigma_mz},
            {Instrument::QUAD, 1, sigma_mz});
        Grid::splash(mesh, 0.0, 0.0, 1, sigma_rt);
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
        for (unsigned int j = 0; j < 7; j++) {
            for (unsigned int i = 0; i < 7; ++i) {
                CHECK(expected_weights[i + 7 * j] ==
                      std::to_string(mesh.value_at(i, j).value()));
            }
        }
    }

    // Point is totally outside the grid.
    SUBCASE("Splatting outside the grid fails") {
        auto mesh = Mesh(
            {7, 7},
            {-3.0 * sigma_rt, 3.0 * sigma_rt, -3.0 * sigma_mz, 3.0 * sigma_mz},
            {Instrument::QUAD, 1, sigma_mz});
        CHECK_FALSE(Grid::splash(mesh, -20.0, -20.0, 1, sigma_rt));
        CHECK_FALSE(Grid::splash(mesh, 20.0, 20.0, 1, sigma_rt));
        CHECK_FALSE(Grid::splash(mesh, -20.0, 20.0, 1, sigma_rt));
        CHECK_FALSE(Grid::splash(mesh, 20.0, +20.0, 1, sigma_rt));
    }

    SUBCASE("Splat outside but range falls inside") {
        auto mesh = Mesh(
            {7, 7},
            {-3.0 * sigma_rt, 3.0 * sigma_rt, -3.0 * sigma_mz, 3.0 * sigma_mz},
            {Instrument::QUAD, 1, sigma_mz});
        CHECK(Grid::splash(mesh, -5.0, 0.0, 1, sigma_rt));
        CHECK(Grid::splash(mesh, 5.0, 0.0, 1, sigma_rt));
        CHECK(Grid::splash(mesh, 0.0, -5.0, 1, sigma_rt));
        CHECK(Grid::splash(mesh, 0.0, 5.0, 1, sigma_rt));
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
        for (unsigned int j = 0; j < 7; j++) {
            for (unsigned int i = 0; i < 7; ++i) {
                CHECK(expected_weights[i + 7 * j] ==
                      std::to_string(mesh.value_at(i, j).value()));
            }
        }
    }
}
