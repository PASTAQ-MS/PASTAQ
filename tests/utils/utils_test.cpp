#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "doctest.h"
#include "mock_stream.hpp"

TEST_CASE("MockStream") {
    SUBCASE("Data access methods <char>(RO)") {
        auto data = std::vector<char>{
            1, 2, 3, 4, 5,  // Row 1
            6, 7, 8, 9, 1,  // Row 2
        };
        MockStream<char> stream(data);
        CHECK(stream.peek() == 1);
        CHECK(stream.good());
        CHECK(stream.tellg() == 0);
        stream.seekg(1, std::ios::cur);
        CHECK(stream.peek() == 2);
        CHECK(stream.good());
        CHECK(stream.tellg() == 1);
        stream.seekg(2, std::ios::cur);
        CHECK(stream.peek() == 4);
        CHECK(stream.good());
        CHECK(stream.tellg() == 3);
        stream.get();
        CHECK(stream.peek() == 5);
        CHECK(stream.good());
        CHECK(stream.tellg() == 4);
        stream.seekg(8, std::ios::beg);
        CHECK(stream.peek() == 9);
        CHECK(stream.good());
        CHECK(stream.tellg() == 8);
        stream.seekg(0, std::ios::beg);
        CHECK(stream.peek() == 1);
        CHECK(stream.good());
        CHECK(stream.tellg() == 0);
        stream.seekg(-1, std::ios::end);
        CHECK(stream.peek() == 1);
        CHECK(stream.good());
        CHECK(stream.tellg() == 9);
    }

    SUBCASE("Data access methods <int>(RO)") {
        auto data = std::vector<int>{
            1, 2, 3, 4, 5,  // Row 1
            6, 7, 8, 9, 1,  // Row 2
        };
        MockStream<int> stream(data);
        // Reading the data element by element.
        int ret = -1;
        for (const auto& e : data) {
            stream.read((char*)&ret, sizeof(int));
            CHECK(ret == e);
            CHECK(stream.good());
        }
        // Attempting to read past the number of elements contained results in
        // stream failure.
        stream.read((char*)&ret, sizeof(int));
        CHECK_FALSE(stream.good());
    }

    SUBCASE("Data access methods <double>(RO)") {
        auto data = std::vector<double>{
            1, 2, 3, 4, 5,  // Row 1
            6, 7, 8, 9, 1,  // Row 2
        };
        MockStream<double> stream(data);
        double ret = -1;
        for (const auto& e : data) {
            stream.read((char*)&ret, sizeof(double));
            CHECK(ret == e);
            CHECK(stream.good());
        }
        // Attempting to read past the number of elements contained results in
        // stream failure.
        stream.read((char*)&ret, sizeof(double));
        CHECK_FALSE(stream.good());
    }

    SUBCASE("Data write methods <int>(RO)") {
        auto data_source = std::vector<int>{
            1, 2, 3, 4, 5,  // Row 1
            6, 7, 8, 9, 1,  // Row 2
        };
        std::vector<int> data_destination(data_source.size());
        MockStream<int> stream(data_destination);
        // Reading the data element by element.
        for (const auto& e : data_source) {
            stream.write((char*)&e, sizeof(int));
        }
        for (size_t i = 0; i < data_source.size(); ++i) {
            CHECK(data_source[i] == data_destination[i]);
        }
        // Attempting to write past the number of elements contained results in
        // stream failure.
        int x = 41;
        stream.write((char*)&x, sizeof(int));
        CHECK_FALSE(stream.good());
    }

    SUBCASE("Data write methods <double>(RO)") {
        auto data_source = std::vector<double>{
            1, 2, 3, 4, 5,  // Row 1
            6, 7, 8, 9, 1,  // Row 2
        };
        std::vector<double> data_destination(data_source.size());
        MockStream<double> stream(data_destination);
        // Reading the data element by element.
        for (const auto& e : data_source) {
            stream.write((char*)&e, sizeof(double));
        }
        for (size_t i = 0; i < data_source.size(); ++i) {
            CHECK(data_source[i] == data_destination[i]);
        }
        // Attempting to write past the number of elements contained results in
        // stream failure.
        double x = 41;
        stream.write((char*)&x, sizeof(double));
        CHECK_FALSE(stream.good());
    }

    SUBCASE("Data read/write methods <int>(RO)") {
        auto data_source = std::vector<int>{
            1, 2, 3, 4, 5,  // Row 1
            6, 7, 8, 9, 1,  // Row 2
        };
        std::vector<int> data_destination(data_source.size());
        MockStream<int> stream(data_destination);
        // Reading the data element by element as well as writing on the
        // destination array. The two pointers for read/write work
        // independently.
        int ret = -1;
        for (const auto& e : data_source) {
            stream.write((char*)&e, sizeof(int));
            stream.read((char*)&ret, sizeof(int));
            CHECK(ret == e);
            CHECK(stream.good());
        }
        for (size_t i = 0; i < data_source.size(); ++i) {
            CHECK(data_source[i] == data_destination[i]);
        }
        // Attempting to write past the number of elements contained results in
        // stream failure.
        int x = 41;
        stream.write((char*)&x, sizeof(int));
        CHECK_FALSE(stream.good());
    }

    SUBCASE("Data read/write methods <double>(RO)") {
        auto data_source = std::vector<double>{
            1, 2, 3, 4, 5,  // Row 1
            6, 7, 8, 9, 1,  // Row 2
        };
        std::vector<double> data_destination(data_source.size());
        MockStream<double> stream(data_destination);
        // Reading the data element by element as well as writing on the
        // destination array. The two pointers for read/write work
        // independently.
        double ret = -1;
        for (const auto& e : data_source) {
            stream.write((char*)&e, sizeof(double));
            stream.read((char*)&ret, sizeof(double));
            CHECK(ret == e);
            CHECK(stream.good());
        }
        for (size_t i = 0; i < data_source.size(); ++i) {
            CHECK(data_source[i] == data_destination[i]);
        }
        // Attempting to write past the number of elements contained results in
        // stream failure.
        double x = 41;
        stream.write((char*)&x, sizeof(double));
        CHECK_FALSE(stream.good());
    }
}
