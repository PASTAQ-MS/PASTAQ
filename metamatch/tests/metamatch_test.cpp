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

std::string file_list =
    R"(file_01.csv 0
file_02.csv 0
file_03.csv 0
file_03.csv 0
file_04.csv 0
file_05.csv 0
file_06.csv 0
file_07.csv 0
file_08.csv 0
file_09.csv 1
file_10.csv 1
file_11.csv 1
file_12.csv 1
file_13.csv 1
file_14.csv 1
file_15.csv 1
file_16.csv 1
file_17.csv 1
)";

// TODO: Add more thourough testing for all functions explored here.
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
    MetaMatch::Parameters parameters = {0.01, 15, 3, 1, 0.6};
    MetaMatch::find_candidates(peaks, parameters);
    MetaMatch::extract_orphans(peaks);
    auto clusters = MetaMatch::reduce_cluster(peaks, 3);
     MetaMatch::write_clusters(std::cout, clusters, parameters);
    //{
        //std::stringstream stream(file_list);
        //auto files = MetaMatch::read_file_list(stream);
        //for (const auto& [file, class_id] : files) {
            //std::cout << "file: " << file << " class_id: " << class_id
                      //<< std::endl;
        //}
    //}

    CHECK(false);
}
