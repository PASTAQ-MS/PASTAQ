#include <iostream>

#include "doctest.h"
#include "tapp-weighted-mesh.hpp"

TEST_CASE("Check if data and weights are being set properly") {
    auto mesh = WeightedMesh({2, 2}, {0.0, 75.0, 200.0, 800.0});
    mesh.set(0,0,10.0,10.0);
    CHECK(mesh.at(0,0) == 100.0);
    CHECK(mesh.weightAt(0,0) == 10.0);
    CHECK(mesh.countsAt(0,0) == 1);
    mesh.set(1,1,1.0,20.0);
    CHECK(mesh.at(1,1) == 20.0);
    CHECK(mesh.weightAt(1,1) == 20.0);
    CHECK(mesh.countsAt(1,1) == 1);
    mesh.set(1,1,2.0,20.0);
    CHECK(mesh.at(1,1) == 60.0);
    CHECK(mesh.weightAt(1,1) == 40.0);
    CHECK(mesh.countsAt(1,1) == 2);
}

TEST_CASE("X/Y indexing") {
    // rt\mz (i,j)  | 0     100   200   300   400   500
    //            0 | (0,0) (1,0) (2,0) (3,0) (4,0) (5,0)
    //           10 | (0,1) (1,1) (2,1) (3,1) (4,1) (5,1)
    //           20 | (0,2) (1,2) (2,2) (3,2) (4,2) (5,2)
    //           30 | (0,3) (1,3) (2,3) (3,3) (4,3) (5,3)
    //           40 | (0,4) (1,4) (2,4) (3,4) (4,4) (5,4)
    auto mesh = WeightedMesh({5, 4}, {0.0, 40.0, 0.0, 500.0});
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
