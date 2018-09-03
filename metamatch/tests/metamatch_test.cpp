#include <iostream>
#include <map>
#include <sstream>

#include "doctest.h"

#include "metamatch/metamatch.hpp"

TEST_CASE("Simple clustering of small peaks") {
    SUBCASE("Peak list is read properly") {
        std::map<std::string, std::string> peak_files = {
            {"file_01.csv",
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
)"},
            {"file_02.csv",
             R"(N X Y Height Volume VCentroid XSigma YSigma Count LocalBkgnd SNVolume SNHeight SNCentroid
0 200.0 100.0 1.0e+08 5.0e+08 4.0e+08 0.01 10 1 1.0e+05 1000 5000 1000
1 300.0 100.0 0.8e+08 3.0e+08 4.0e+08 0.01 10 1 0.8e+05 1000 5000 1000
2 400.0 100.0 0.6e+08 1.0e+08 4.0e+08 0.01 10 1 0.6e+05 1000 5000 1000
3 450.0 200.0 0.8e+08 5.0e+08 4.0e+08 0.01 10 1 0.8e+05 1000 5000 1000
4 650.0 200.0 0.6e+08 5.0e+08 4.0e+08 0.01 10 1 0.6e+05 1000 5000 1000
5 500.0 300.0 0.8e+08 5.0e+08 4.0e+08 0.01 10 1 0.6e+05 1000 5000 1000
)"},
            {"file_03.csv",
             R"(N X Y Height Volume VCentroid XSigma YSigma Count LocalBkgnd SNVolume SNHeight SNCentroid
0 200.0 100.0 1.0e+08 5.0e+08 4.0e+08 0.01 10 1 1.0e+05 1000 5000 1000
1 300.0 100.0 0.8e+08 3.0e+08 4.0e+08 0.01 10 1 0.8e+05 1000 5000 1000
2 400.0 100.0 0.6e+08 1.0e+08 4.0e+08 0.01 10 1 0.6e+05 1000 5000 1000
3 250.0 200.0 1.0e+08 5.0e+08 4.0e+08 0.01 10 1 1.0e+05 1000 5000 1000
4 450.0 200.0 0.8e+08 5.0e+08 4.0e+08 0.01 10 1 0.8e+05 1000 5000 1000
5 650.0 200.0 0.6e+08 5.0e+08 4.0e+08 0.01 10 1 0.6e+05 1000 5000 1000
6 300.0 300.0 0.8e+08 5.0e+08 4.0e+08 0.01 10 1 0.6e+05 1000 5000 1000
)"},
        };

        std::string file_list =
            R"(file_01.csv 0
file_02.csv 0
file_03.csv 0
file_01.csv 1
file_02.csv 1
file_03.csv 1
)";

        // Read file list.
        std::stringstream stream(file_list);
        auto files = MetaMatch::read_file_list(stream);
        CHECK(files.size() == 6);
        if (files.size() == 6) {
            CHECK(files[0].first == "file_01.csv");
            CHECK(files[0].second == 0);
            CHECK(files[1].first == "file_02.csv");
            CHECK(files[1].second == 0);
            CHECK(files[2].first == "file_03.csv");
            CHECK(files[2].second == 0);
            CHECK(files[3].first == "file_01.csv");
            CHECK(files[3].second == 1);
            CHECK(files[4].first == "file_02.csv");
            CHECK(files[4].second == 1);
            CHECK(files[5].first == "file_03.csv");
            CHECK(files[5].second == 1);
            // Read peaks into MetaMatch::Peak array.
            int i = 0;
            std::vector<MetaMatch::Peak> peaks;
            for (const auto& [file, class_id] : files) {
                std::vector<Centroid::Peak> file_peaks;
                std::cout << "file: " << file << " class_id: " << class_id
                          << std::endl;
                std::stringstream stream(peak_files[file]);
                Centroid::Files::Csv::read_peaks(stream, &file_peaks);
                for (const auto& peak : file_peaks) {
                    peaks.push_back({peak, i, class_id, -1, peak.mz, peak.rt});
                }
                ++i;
            }
            CHECK(peaks.size() == 44);
            if (peaks.size() == 44) {
                MetaMatch::Parameters parameters = {0.01, 15, 0.6};
                MetaMatch::find_candidates(peaks, parameters);
                MetaMatch::extract_orphans(peaks);
                CHECK(peaks.size() == 34);
                auto clusters = MetaMatch::reduce_cluster(peaks, files.size());
                CHECK(clusters.size() == 5);
            }
        }
    }
    // CHECK(false);
}
