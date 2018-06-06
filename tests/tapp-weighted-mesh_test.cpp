#include <iostream>

#include "doctest.h"
#include "tapp-weighted-mesh.hpp"

TEST_CASE("Bounds check on WeightedMesh when using WeightedMesh::at") {
    SUBCASE("An empty mesh should always return std::nullopt") {
        auto mesh = WeightedMesh();
        CHECK(mesh.at(0, 0) == std::nullopt);
        CHECK(mesh.at(0, 1) == std::nullopt);
        CHECK(mesh.at(1, 0) == std::nullopt);
        CHECK(mesh.at(1, 1) == std::nullopt);
    };
    SUBCASE(
        "An mesh should return it's value if inside bounds or std::nullopt") {
        auto mesh = WeightedMesh({4, 4}, {0.0, 60.0, 80.5, 1000.0});
        CHECK(mesh.at(0, 0) == 0.0);
        CHECK(mesh.at(0, 0) == 0.0);
        CHECK(mesh.at(0, 1) == 0.0);
        CHECK(mesh.at(1, 0) == 0.0);
        CHECK(mesh.at(1, 1) == 0.0);
        CHECK(mesh.at(3, 3) == 0.0);
        CHECK(mesh.at(3, 4) == std::nullopt);
        CHECK(mesh.at(4, 3) == std::nullopt);
    };
}

TEST_CASE(
    "Real world dimensions fetching via WeightedMesh::mz and "
    "WeightedMesh::rt") {
    SUBCASE("An empty mesh should always return std::nullopt") {
        auto mesh = WeightedMesh();
        CHECK(mesh.mz(0) == std::nullopt);
        CHECK(mesh.mz(1) == std::nullopt);
        CHECK(mesh.mz(2) == std::nullopt);
        CHECK(mesh.rt(0) == std::nullopt);
        CHECK(mesh.rt(1) == std::nullopt);
        CHECK(mesh.rt(2) == std::nullopt);
    };
    SUBCASE("Proper interpolation if inside bounds or std::nullopt") {
        auto mesh = WeightedMesh({4, 4}, {0.0, 75.0, 200.0, 800.0});
        CHECK(mesh.mz(0) == 200.0);
        CHECK(mesh.mz(1) == 400.0);
        CHECK(mesh.mz(2) == 600.0);
        CHECK(mesh.mz(3) == 800.0);
        CHECK(mesh.mz(4) == std::nullopt);
        CHECK(mesh.rt(0) == 0.0);
        CHECK(mesh.rt(1) == 25.0);
        CHECK(mesh.rt(2) == 50.0);
        CHECK(mesh.rt(3) == 75.0);
        CHECK(mesh.rt(4) == std::nullopt);
    };
}

TEST_CASE("Getters for mDimensions and mBounds") {
    SUBCASE("With an empty mesh") {
        auto mesh = WeightedMesh();
        auto dimensions = mesh.dim();
        auto bounds = mesh.bounds();
        CHECK(dimensions.n == 0);
        CHECK(dimensions.m == 0);
        CHECK(bounds.minRt == 0.0);
        CHECK(bounds.maxRt == 0.0);
        CHECK(bounds.minMz == 0.0);
        CHECK(bounds.maxMz == 0.0);
    };
    SUBCASE("With a non empty mesh") {
        auto mesh = WeightedMesh({4, 10}, {0.0, 75.0, 200.0, 800.0});
        auto dimensions = mesh.dim();
        auto bounds = mesh.bounds();
        CHECK(dimensions.n == 4);
        CHECK(dimensions.m == 10);
        CHECK(bounds.minRt == 0.0);
        CHECK(bounds.maxRt == 75.0);
        CHECK(bounds.minMz == 200.0);
        CHECK(bounds.maxMz == 800.0);
    };
}

TEST_CASE("Check if data is being set properly") {
    auto mesh = WeightedMesh({2, 2}, {0.0, 75.0, 200.0, 800.0});
    mesh.printAll();
    mesh.set(0,0,10,10);
    mesh.set(1,1,1,20);
    mesh.set(1,1,2,20);
    mesh.printAll();
    CHECK(1 == 0);
}
