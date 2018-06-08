#include <iostream>
#include <string>

#include "doctest.h"
#include "tapp-weighted-mesh.hpp"

TEST_CASE("X/Y indexing") {
    // rt\mz (i,j)  | 0     100   200   300   400   500
    //            0 | (0,0) (1,0) (2,0) (3,0) (4,0) (5,0)
    //           10 | (0,1) (1,1) (2,1) (3,1) (4,1) (5,1)
    //           20 | (0,2) (1,2) (2,2) (3,2) (4,2) (5,2)
    //           30 | (0,3) (1,3) (2,3) (3,3) (4,3) (5,3)
    //           40 | (0,4) (1,4) (2,4) (3,4) (4,4) (5,4)
    auto mesh = WeightedMesh({6, 5}, {0.0, 40.0, 0.0, 500.0});
    CHECK(mesh.x_index(0.0) == 0);
    CHECK(mesh.x_index(10.0) == 0);
    CHECK(mesh.x_index(40.0) == 0);
    CHECK(mesh.x_index(60.0) == 0);
    CHECK(mesh.x_index(80.0) == 0);
    CHECK(mesh.x_index(100.0) == 1);
    CHECK(mesh.x_index(250.0) == 2);
    CHECK(mesh.x_index(333.0) == 3);
    CHECK(mesh.x_index(499.0) == 4);
    CHECK(mesh.x_index(500.0) == 5);
    CHECK(mesh.x_index(501.0) == 5);
    CHECK(mesh.x_index(599.0) == 5);
    CHECK(mesh.x_index(600.0) == std::nullopt);
    CHECK(mesh.x_index(601.0) == std::nullopt);
    CHECK(mesh.y_index(0.0) == 0);
    CHECK(mesh.y_index(1.0) == 0);
    CHECK(mesh.y_index(4.0) == 0);
    CHECK(mesh.y_index(6.0) == 0);
    CHECK(mesh.y_index(8.0) == 0);
    CHECK(mesh.y_index(10.0) == 1);
    CHECK(mesh.y_index(25.0) == 2);
    CHECK(mesh.y_index(33.0) == 3);
    CHECK(mesh.y_index(40.0) == 4);
    CHECK(mesh.y_index(49.0) == 4);
    CHECK(mesh.y_index(50.0) == std::nullopt);
    CHECK(mesh.y_index(59.0) == std::nullopt);
}

TEST_CASE("Gaussian splatting") {
    double sigma_mz = 1.0;
    double sigma_rt = 1.0;

    // For testing purposes we will set the bounds to [-3 * sigma, +3 * sigma]
    // on both dimensions, even if there is no such thing as negative retention
    // time or mass.
    auto mesh = WeightedMesh({7, 7}, {-3.0 * sigma_rt, 3.0 * sigma_rt,
                                      -3.0 * sigma_mz, 3.0 * sigma_mz});
    mesh.splash(0.0, 0.0, 10, sigma_mz, sigma_rt);
    CHECK(mesh.at(3, 3) == 10);
    CHECK(mesh.weight_at(3, 3) == 1);
    CHECK(std::to_string(mesh.weight_at(2, 3).value()) == "0.606531");
    CHECK(std::to_string(mesh.weight_at(3, 2).value()) == "0.606531");
    CHECK(std::to_string(mesh.weight_at(3, 4).value()) == "0.606531");
    CHECK(std::to_string(mesh.weight_at(4, 3).value()) == "0.606531");
    CHECK(std::to_string(mesh.weight_at(1, 3).value()) == "0.135335");
    CHECK(std::to_string(mesh.weight_at(3, 1).value()) == "0.135335");
    CHECK(std::to_string(mesh.weight_at(3, 5).value()) == "0.135335");
    CHECK(std::to_string(mesh.weight_at(5, 3).value()) == "0.135335");
    CHECK(std::to_string(mesh.at(2, 3).value()) == "6.065307");
    CHECK(std::to_string(mesh.at(3, 2).value()) == "6.065307");
    CHECK(std::to_string(mesh.at(3, 4).value()) == "6.065307");
    CHECK(std::to_string(mesh.at(4, 3).value()) == "6.065307");
    CHECK(std::to_string(mesh.at(1, 3).value()) == "1.353353");
    CHECK(std::to_string(mesh.at(3, 1).value()) == "1.353353");
    CHECK(std::to_string(mesh.at(3, 5).value()) == "1.353353");
    CHECK(std::to_string(mesh.at(5, 3).value()) == "1.353353");
}
