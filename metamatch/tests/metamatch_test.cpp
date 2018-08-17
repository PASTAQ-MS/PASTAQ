#include <iostream>
#include <sstream>

#include "doctest.h"

#include "metamatch/metamatch.hpp"

const char* file_01 =
    R"(N X Y Height Volume VCentroid XSigma YSigma Count LocalBkgnd SNVolume SNHeight SNCentroid
0 200.0 100.0 1.0e+08 5.0e+08 4.0e+08 0.01 10 1 1.0e+05 1000 5000 1000
1 300.0 100.0 0.8e+08 3.0e+08 4.0e+08 0.01 10 1 0.8e+05 1000 5000 1000
2 400.0 100.0 0.6e+08 1.0e+08 4.0e+08 0.01 10 1 0.6e+05 1000 5000 1000
3 250.0 200.0 1.0e+08 5.0e+08 4.0e+08 0.01 10 1 1.0e+05 1000 5000 1000
4 450.0 200.0 0.8e+08 5.0e+08 4.0e+08 0.01 10 1 0.8e+05 1000 5000 1000
5 500.0 200.0 1.0e+08 5.0e+08 4.0e+08 0.01 10 1 1.0e+05 1000 5000 1000
6 650.0 200.0 0.6e+08 5.0e+08 4.0e+08 0.01 10 1 0.6e+05 1000 5000 1000
7 800.0 300.0 0.8e+08 5.0e+08 4.0e+08 0.01 10 1 0.6e+05 1000 5000 1000
8 200.0 300.0 0.5e+08 5.0e+08 4.0e+08 0.01 10 1 0.6e+05 1000 5000 1000
)";

TEST_CASE("DEBUG") {
    std::stringstream stream(file_01);
    std::string line;
    int i = 0;
    while (std::getline(stream, line)) {
        std::cout << "line [" << i << "]: " << line << std::endl;
        ++i;
    }
    //CHECK(false);
}
