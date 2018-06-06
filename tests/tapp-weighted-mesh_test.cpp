#include <iostream>

#include "doctest.h"
#include "tapp-weighted-mesh.hpp"

TEST_CASE("Check if data and weights are being set properly") {
    auto mesh = WeightedMesh({2, 2}, {0.0, 75.0, 200.0, 800.0});
    mesh.set(0,0,10,10);
    CHECK(mesh.at(0,0) == 100);
    CHECK(mesh.weightAt(0,0) == 10);
    CHECK(mesh.countsAt(0,0) == 1);
    mesh.set(1,1,1,20);
    CHECK(mesh.at(1,1) == 20);
    CHECK(mesh.weightAt(1,1) == 20);
    CHECK(mesh.countsAt(1,1) == 1);
    mesh.set(1,1,2,20);
    CHECK(mesh.at(1,1) == 60);
    CHECK(mesh.weightAt(1,1) == 40);
    CHECK(mesh.countsAt(1,1) == 2);
}
