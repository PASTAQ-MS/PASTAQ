#include <iostream>

#include "doctest.h"
#include "tapp-mesh.hpp"

TEST_CASE("Bounds check on Mesh when using Mesh::at") {
    SUBCASE("An empty mesh should always return std::nullopt") {
        auto mesh = Mesh();
        CHECK(mesh.value_at(0, 0) == std::nullopt);
        CHECK(mesh.value_at(0, 1) == std::nullopt);
        CHECK(mesh.value_at(1, 0) == std::nullopt);
        CHECK(mesh.value_at(1, 1) == std::nullopt);
    };
    SUBCASE(
        "An mesh should return it's value if inside bounds or std::nullopt") {
        auto mesh = Mesh({4, 4}, {0.0, 60.0, 80.5, 1000.0});
        CHECK(mesh.value_at(0, 0) == 0.0);
        CHECK(mesh.value_at(0, 0) == 0.0);
        CHECK(mesh.value_at(0, 1) == 0.0);
        CHECK(mesh.value_at(1, 0) == 0.0);
        CHECK(mesh.value_at(1, 1) == 0.0);
        CHECK(mesh.value_at(3, 3) == 0.0);
        CHECK(mesh.value_at(3, 4) == std::nullopt);
        CHECK(mesh.value_at(4, 3) == std::nullopt);
    };
}

TEST_CASE("Setting values at given indices should work properly") {
    auto mesh = Mesh({2, 2}, {0.0, 75.0, 200.0, 800.0});
    mesh.set_value(0, 0, 10.0);
    CHECK(mesh.value_at(0, 0) == 10.0);
    mesh.set_value(1, 1, 1.0);
    CHECK(mesh.value_at(1, 1) == 1.0);
    mesh.set_value(1, 1, 2.0);
    CHECK(mesh.value_at(1, 1) == 2.0);
    mesh.set_value(0, 1, 20.0);
    CHECK(mesh.value_at(0, 1) == 20.0);
}

TEST_CASE("Real world dimensions fetching via Mesh::mz and Mesh::rt") {
    SUBCASE("An empty mesh should always return std::nullopt") {
        auto mesh = Mesh();
        CHECK(mesh.mz_at(0) == std::nullopt);
        CHECK(mesh.mz_at(1) == std::nullopt);
        CHECK(mesh.mz_at(2) == std::nullopt);
        CHECK(mesh.rt_at(0) == std::nullopt);
        CHECK(mesh.rt_at(1) == std::nullopt);
        CHECK(mesh.rt_at(2) == std::nullopt);
    };
    SUBCASE("Proper interpolation if inside bounds or std::nullopt") {
        auto mesh = Mesh({4, 4}, {0.0, 75.0, 200.0, 800.0});
        CHECK(mesh.mz_at(0) == 200.0);
        CHECK(mesh.mz_at(1) == 400.0);
        CHECK(mesh.mz_at(2) == 600.0);
        CHECK(mesh.mz_at(3) == 800.0);
        CHECK(mesh.mz_at(4) == std::nullopt);
        CHECK(mesh.rt_at(0) == 0.0);
        CHECK(mesh.rt_at(1) == 25.0);
        CHECK(mesh.rt_at(2) == 50.0);
        CHECK(mesh.rt_at(3) == 75.0);
        CHECK(mesh.rt_at(4) == std::nullopt);
    };
}

TEST_CASE("Getters for mDimensions and mBounds") {
    SUBCASE("With an empty mesh") {
        auto mesh = Mesh();
        auto dimensions = mesh.dim();
        auto bounds = mesh.bounds();
        CHECK(dimensions.n == 0);
        CHECK(dimensions.m == 0);
        CHECK(bounds.min_rt == 0.0);
        CHECK(bounds.max_rt == 0.0);
        CHECK(bounds.min_mz == 0.0);
        CHECK(bounds.max_mz == 0.0);
    };
    SUBCASE("With a non empty mesh") {
        auto mesh = Mesh({4, 10}, {0.0, 75.0, 200.0, 800.0});
        auto dimensions = mesh.dim();
        auto bounds = mesh.bounds();
        CHECK(dimensions.n == 4);
        CHECK(dimensions.m == 10);
        CHECK(bounds.min_rt == 0.0);
        CHECK(bounds.max_rt == 75.0);
        CHECK(bounds.min_mz == 200.0);
        CHECK(bounds.max_mz == 800.0);
    };
}
