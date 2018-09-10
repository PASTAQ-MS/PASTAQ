#include <map>
#include <sstream>

#include "doctest.h"

#include "metamatch/metamatch.hpp"

TEST_CASE("Clustering small peak list") {
    std::map<std::string, std::string> peak_files = {
        {
            "file_01.csv",
            "N X Y Height Volume VCentroid XSigma YSigma Count LocalBkgnd "
            "SNVolume SNHeight SNCentroid\n"
            "0 200.0 100.0 1 1 1 0.01 10 1 1 1000 5000 1000\n"
            "1 300.0 100.0 1 1 1 0.01 10 1 1 1000 5000 1000\n"
            "2 400.0 100.0 1 1 1 0.01 10 1 1 1000 5000 1000\n"
            "3 250.0 200.0 1 1 1 0.01 10 1 1 1000 5000 1000\n"
            "4 450.0 200.0 1 1 1 0.01 10 1 1 1000 5000 1000\n"
            "5 500.0 200.0 1 1 1 0.01 10 1 1 1000 5000 1000\n"
            "6 650.0 200.0 1 1 1 0.01 10 1 1 1000 5000 1000\n"
            "7 800.0 300.0 1 1 1 0.01 10 1 1 1000 5000 1000\n"
            "8 200.0 300.0 1 1 1 0.01 10 1 1 1000 5000 1000\n",
        },
        {
            "file_02.csv",
            "N X Y Height Volume VCentroid XSigma YSigma Count LocalBkgnd "
            "SNVolume SNHeight SNCentroid\n"
            "0 200.0 100.0 1 1 1 0.01 10 1 1 1000 5000 1000\n"
            "1 300.0 100.0 1 1 1 0.01 10 1 1 1000 5000 1000\n"
            "2 400.0 100.0 1 1 1 0.01 10 1 1 1000 5000 1000\n"
            "3 450.0 200.0 1 1 1 0.01 10 1 1 1000 5000 1000\n"
            "4 650.0 200.0 1 1 1 0.01 10 1 1 1000 5000 1000\n"
            "5 500.0 300.0 1 1 1 0.01 10 1 1 1000 5000 1000\n",
        },
        {
            "file_03.csv",
            "N X Y Height Volume VCentroid XSigma YSigma Count LocalBkgnd "
            "SNVolume SNHeight SNCentroid\n"
            "0 200.0 100.0 1 1 1 0.01 10 1 1 1000 5000 1000\n"
            "1 300.0 100.0 1 1 1 0.01 10 1 1 1000 5000 1000\n"
            "2 400.0 100.0 1 1 1 0.01 10 1 1 1000 5000 1000\n"
            "3 250.0 200.0 1 1 1 0.01 10 1 1 1000 5000 1000\n"
            "4 450.0 200.0 1 1 1 0.01 10 1 1 1000 5000 1000\n"
            "5 650.0 200.0 1 1 1 0.01 10 1 1 1000 5000 1000\n"
            "6 300.0 300.0 1 1 1 0.01 10 1 1 1000 5000 1000\n",
        },
    };

    SUBCASE("Peak list 01") {
        std::vector<std::pair<std::string, size_t>> files = {
            {"file_01.csv", 0},
            {"file_02.csv", 0},
            {"file_03.csv", 0},
        };
        // Read peaks into MetaMatch::Peak array.
        size_t i = 0;
        std::vector<MetaMatch::Peak> peaks;
        for (const auto& [file, class_id] : files) {
            std::vector<Centroid::Peak> file_peaks;
            std::stringstream stream(peak_files[file]);
            Centroid::Files::Csv::read_peaks(stream, &file_peaks);
            for (const auto& peak : file_peaks) {
                peaks.push_back({peak, i, class_id, -1, peak.mz, peak.rt});
            }
            ++i;
        }
        CHECK(peaks.size() == 22);
        if (peaks.size() == 22) {
            MetaMatch::Parameters parameters = {0.01, 15, {{0, 3, 2}}};
            MetaMatch::find_clusters(peaks, parameters);
            MetaMatch::extract_orphans(peaks);
            CHECK(peaks.size() == 17);
            std::vector<MetaMatch::Cluster> expected_clusters = {
                {0, 200, 100, {}}, {1, 250, 200, {}}, {2, 300, 100, {}},
                {3, 400, 100, {}}, {4, 450, 200, {}},
            };
            auto clusters = MetaMatch::reduce_cluster(peaks, files.size());
            CHECK(clusters.size() == expected_clusters.size());
            if (clusters.size() == expected_clusters.size()) {
                for (size_t i = 0; i < clusters.size(); ++i) {
                    CHECK(clusters[i].id == expected_clusters[i].id);
                    CHECK(clusters[i].mz == expected_clusters[i].mz);
                    CHECK(clusters[i].rt == expected_clusters[i].rt);
                }
            }
        }
    }

    SUBCASE("Peak list 02") {
        std::vector<std::pair<std::string, size_t>> files = {
            {"file_01.csv", 0}, {"file_02.csv", 0}, {"file_03.csv", 0},
            {"file_01.csv", 1}, {"file_02.csv", 1}, {"file_03.csv", 1},
        };
        // Read peaks into MetaMatch::Peak array.
        size_t i = 0;
        std::vector<MetaMatch::Peak> peaks;
        for (const auto& [file, class_id] : files) {
            std::vector<Centroid::Peak> file_peaks;
            std::stringstream stream(peak_files[file]);
            Centroid::Files::Csv::read_peaks(stream, &file_peaks);
            for (const auto& peak : file_peaks) {
                peaks.push_back({peak, i, class_id, -1, peak.mz, peak.rt});
            }
            ++i;
        }
        CHECK(peaks.size() == 44);
        if (peaks.size() == 44) {
            MetaMatch::Parameters parameters = {
                0.01, 15, {{0, 3, 2}, {1, 3, 2}}};
            MetaMatch::find_clusters(peaks, parameters);
            MetaMatch::extract_orphans(peaks);
            CHECK(peaks.size() == 34);
            std::vector<MetaMatch::Cluster> expected_clusters = {
                {0, 200, 100, {}}, {1, 250, 200, {}}, {2, 300, 100, {}},
                {3, 400, 100, {}}, {4, 450, 200, {}},
            };
            auto clusters = MetaMatch::reduce_cluster(peaks, files.size());
            CHECK(clusters.size() == expected_clusters.size());
            if (clusters.size() == expected_clusters.size()) {
                for (size_t i = 0; i < clusters.size(); ++i) {
                    CHECK(clusters[i].id == expected_clusters[i].id);
                    CHECK(clusters[i].mz == expected_clusters[i].mz);
                    CHECK(clusters[i].rt == expected_clusters[i].rt);
                }
            }
        }
    }
}
