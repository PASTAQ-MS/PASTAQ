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

const char* file_02 =
    R"(N X Y Height Volume VCentroid XSigma YSigma Count LocalBkgnd SNVolume SNHeight SNCentroid
0 200.0 100.0 1.0e+08 5.0e+08 4.0e+08 0.01 10 1 1.0e+05 1000 5000 1000
1 300.0 100.0 0.8e+08 3.0e+08 4.0e+08 0.01 10 1 0.8e+05 1000 5000 1000
2 400.0 100.0 0.6e+08 1.0e+08 4.0e+08 0.01 10 1 0.6e+05 1000 5000 1000
3 450.0 200.0 0.8e+08 5.0e+08 4.0e+08 0.01 10 1 0.8e+05 1000 5000 1000
4 650.0 200.0 0.6e+08 5.0e+08 4.0e+08 0.01 10 1 0.6e+05 1000 5000 1000
5 500.0 300.0 0.8e+08 5.0e+08 4.0e+08 0.01 10 1 0.6e+05 1000 5000 1000
)";

const char* file_03 =
    R"(N X Y Height Volume VCentroid XSigma YSigma Count LocalBkgnd SNVolume SNHeight SNCentroid
0 200.0 100.0 1.0e+08 5.0e+08 4.0e+08 0.01 10 1 1.0e+05 1000 5000 1000
1 300.0 100.0 0.8e+08 3.0e+08 4.0e+08 0.01 10 1 0.8e+05 1000 5000 1000
2 400.0 100.0 0.6e+08 1.0e+08 4.0e+08 0.01 10 1 0.6e+05 1000 5000 1000
3 250.0 200.0 1.0e+08 5.0e+08 4.0e+08 0.01 10 1 1.0e+05 1000 5000 1000
4 450.0 200.0 0.8e+08 5.0e+08 4.0e+08 0.01 10 1 0.8e+05 1000 5000 1000
5 650.0 200.0 0.6e+08 5.0e+08 4.0e+08 0.01 10 1 0.6e+05 1000 5000 1000
6 300.0 300.0 0.8e+08 5.0e+08 4.0e+08 0.01 10 1 0.6e+05 1000 5000 1000
)";

TEST_CASE("DEBUG") {
    std::vector<const char*> files = {file_01, file_02, file_03};
    // Read the files.
    std::vector<std::vector<Centroid::Peak>> file_peaks(files.size());
    std::vector<MetaMatch::Peak> peaks;
    for (size_t i = 0; i < files.size(); ++i) {
        auto file = files[i];
        std::stringstream stream(file);
        std::string line;
        Centroid::Files::Csv::read_peaks(stream, &file_peaks[i]);
        for (const auto& peak : file_peaks[i]) {
            peaks.push_back({peak, i, 0, -1, peak.mz, peak.rt});
        }
    }

    // ...
    MetaMatch::find_candidates(peaks, 0.01, 15, 3, 1, 0.6);

    CHECK(false);
}
