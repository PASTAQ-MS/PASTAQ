#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <streambuf>

#include "doctest.h"
#include "tapp-mesh.hpp"

TEST_CASE("Bounds check on RegularMesh when using RegularMesh::at") {
    SUBCASE("An empty mesh should always return std::nullopt") {
        auto mesh = RegularMesh({0, 0}, {0.0, 60.0, 80.5, 1000.0},
                                Instrument::QUAD, {200, 0.01, 1.0});
        CHECK(mesh.value_at(0, 0) == std::nullopt);
        CHECK(mesh.value_at(0, 1) == std::nullopt);
        CHECK(mesh.value_at(1, 0) == std::nullopt);
        CHECK(mesh.value_at(1, 1) == std::nullopt);
    };
    SUBCASE(
        "An mesh should return it's value if inside bounds or std::nullopt") {
        auto mesh = RegularMesh({4, 4}, {0.0, 60.0, 80.5, 1000.0},
                                Instrument::QUAD, {200, 0.01, 1.0});
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
    auto mesh = RegularMesh({2, 2}, {0.0, 75.0, 200.0, 800.0}, Instrument::QUAD,
                            {200, 0.01, 1.0});
    mesh.set_value(0, 0, 10.0);
    CHECK(mesh.value_at(0, 0) == 10.0);
    mesh.set_value(1, 1, 1.0);
    CHECK(mesh.value_at(1, 1) == 1.0);
    mesh.set_value(1, 1, 2.0);
    CHECK(mesh.value_at(1, 1) == 2.0);
    mesh.set_value(0, 1, 20.0);
    CHECK(mesh.value_at(0, 1) == 20.0);
}

TEST_CASE(
    "Real world dimensions fetching via RegularMesh::mz and RegularMesh::rt") {
    SUBCASE("An empty mesh should always return std::nullopt") {
        auto mesh = RegularMesh({0, 0}, {0.0, 75.0, 200.0, 800.0},
                                Instrument::QUAD, {200, 0.01, 1.0});
        CHECK(mesh.mz_at(0) == std::nullopt);
        CHECK(mesh.mz_at(1) == std::nullopt);
        CHECK(mesh.mz_at(2) == std::nullopt);
        CHECK(mesh.rt_at(0) == std::nullopt);
        CHECK(mesh.rt_at(1) == std::nullopt);
        CHECK(mesh.rt_at(2) == std::nullopt);
    };
    SUBCASE("Proper interpolation if inside bounds or std::nullopt") {
        auto mesh = RegularMesh({4, 4}, {0.0, 75.0, 200.0, 800.0},
                                Instrument::QUAD, {200, 0.01, 1.0});
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
        auto mesh = RegularMesh({0, 0}, {0.0, 0.0, 0.0, 0.0}, Instrument::QUAD,
                                {200, 0.01, 1.0});
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
        auto mesh = RegularMesh({4, 10}, {0.0, 75.0, 200.0, 800.0},
                                Instrument::QUAD, {200, 0.01, 1.0});
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

TEST_CASE("Loading data from the dat file") {
    // auto mesh = RegularMesh({}, {}, Instrument::QUAD, {});
    // auto data = std::vector<double>{
    // 1.0, 2.0, 3.0, 4.0, 5.0, // Row 1
    // 6.0, 7.0, 8.0, 9.0, 0.0, // Row 2
    //};
    auto mesh = RegularMesh({4, 10}, {0.0, 75.0, 200.0, 800.0},
                            Instrument::QUAD, {200, 0.01, 1.0});
    mesh.set_value(0, 0, 10.0);
    mesh.set_value(1, 1, 1.0);
    mesh.set_value(1, 1, 2.0);
    mesh.set_value(0, 1, 20.0);
    std::ofstream fileout("foo.dat", std::ios::out | std::ios::binary);
    CHECK(mesh.write_dat(fileout));
    // fileout << 1.0;

    // std::istringstream stream;
    // std::iostream *sbuff;
    // std::copy(data.begin(), data.end(), std::istream_iterator<double>(sbuff)
    // );
}
