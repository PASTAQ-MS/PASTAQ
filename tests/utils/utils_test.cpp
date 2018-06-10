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
            std::cout << stream.eof() << std::endl;
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
            std::cout << stream.eof() << std::endl;
        }
        // Attempting to read past the number of elements contained results in
        // stream failure.
        stream.read((char*)&ret, sizeof(double));
        CHECK_FALSE(stream.good());
    }
}
