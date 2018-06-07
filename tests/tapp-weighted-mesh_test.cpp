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
    CHECK(mesh.xIndex(0.0) == 0);
    CHECK(mesh.xIndex(10.0) == 0);
    CHECK(mesh.xIndex(40.0) == 0);
    CHECK(mesh.xIndex(60.0) == 0);
    CHECK(mesh.xIndex(80.0) == 0);
    CHECK(mesh.xIndex(100.0) == 1);
    CHECK(mesh.xIndex(250.0) == 2);
    CHECK(mesh.xIndex(333.0) == 3);
    CHECK(mesh.xIndex(499.0) == 4);
    CHECK(mesh.xIndex(500.0) == 5);
    CHECK(mesh.xIndex(501.0) == 5);
    CHECK(mesh.xIndex(599.0) == 5);
    CHECK(mesh.xIndex(600.0) == std::nullopt);
    CHECK(mesh.xIndex(601.0) == std::nullopt);
    CHECK(mesh.yIndex(0.0) == 0);
    CHECK(mesh.yIndex(1.0) == 0);
    CHECK(mesh.yIndex(4.0) == 0);
    CHECK(mesh.yIndex(6.0) == 0);
    CHECK(mesh.yIndex(8.0) == 0);
    CHECK(mesh.yIndex(10.0) == 1);
    CHECK(mesh.yIndex(25.0) == 2);
    CHECK(mesh.yIndex(33.0) == 3);
    CHECK(mesh.yIndex(40.0) == 4);
    CHECK(mesh.yIndex(49.0) == 4);
    CHECK(mesh.yIndex(50.0) == std::nullopt);
    CHECK(mesh.yIndex(59.0) == std::nullopt);
}

TEST_CASE("Gaussian splatting") {
    double sigmaMz = 1.0;
    double sigmaRt = 1.0;

    // For testing purposes we will set the bounds to [-3 * sigma, +3 * sigma]
    // on both dimensions, even if there is no such thing as negative retention
    // time or mass.
    auto mesh = WeightedMesh(
        {7, 7}, {-3.0 * sigmaRt, 3.0 * sigmaRt, -3.0 * sigmaMz, 3.0 * sigmaMz});
    mesh.splash(0.0, 0.0, 10, sigmaMz, sigmaRt);
    CHECK(mesh.at(3, 3) == 10);
    CHECK(mesh.weightAt(3, 3) == 1);
    CHECK(std::to_string(mesh.weightAt(2, 3).value()) == "0.606531");
    CHECK(std::to_string(mesh.weightAt(3, 2).value()) == "0.606531");
    CHECK(std::to_string(mesh.weightAt(3, 4).value()) == "0.606531");
    CHECK(std::to_string(mesh.weightAt(4, 3).value()) == "0.606531");
    CHECK(std::to_string(mesh.weightAt(1, 3).value()) == "0.135335");
    CHECK(std::to_string(mesh.weightAt(3, 1).value()) == "0.135335");
    CHECK(std::to_string(mesh.weightAt(3, 5).value()) == "0.135335");
    CHECK(std::to_string(mesh.weightAt(5, 3).value()) == "0.135335");
    CHECK(std::to_string(mesh.at(2, 3).value()) == "6.065307");
    CHECK(std::to_string(mesh.at(3, 2).value()) == "6.065307");
    CHECK(std::to_string(mesh.at(3, 4).value()) == "6.065307");
    CHECK(std::to_string(mesh.at(4, 3).value()) == "6.065307");
    CHECK(std::to_string(mesh.at(1, 3).value()) == "1.353353");
    CHECK(std::to_string(mesh.at(3, 1).value()) == "1.353353");
    CHECK(std::to_string(mesh.at(3, 5).value()) == "1.353353");
    CHECK(std::to_string(mesh.at(5, 3).value()) == "1.353353");
}
