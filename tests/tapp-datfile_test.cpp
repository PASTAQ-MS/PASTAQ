#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <streambuf>

#include "doctest.h"
#include "tapp-datfile.hpp"

// Helper class to use a vector to mock istreams/ostreams. Note that the
// lifetime of the VectorStream class is not tied to the object we have mapped.
// If the vector get's deallocated or it's dimensions change, it will result in
// an error when using the stream.
struct VectorStream : std::streambuf {
    VectorStream(std::vector<double> &data) {
        auto begin = reinterpret_cast<char *>(&data[0]);
        auto end =
            reinterpret_cast<char *>(&data[0]) + sizeof(double) * data.size();
        setg(begin, begin, end);
        setp(begin, end);
    }
    VectorStream(std::vector<char> &data) {
        auto begin = reinterpret_cast<char *>(&data[0]);
        auto end =
            reinterpret_cast<char *>(&data[0]) + sizeof(char) * data.size();
        setg(begin, begin, end);
        setp(begin, end);
    }
};

TEST_CASE("BLAH") {
    auto data = std::vector<char>{
        1, 2, 3, 4, 5,    // Row 1
        6, 7, 8, 9, 10,  // Row 2
    };

    VectorStream vector_stream(data);
    std::istream stream(&vector_stream);
    DatFile::load_parameters(stream);
    CHECK(1 == 0);
}
