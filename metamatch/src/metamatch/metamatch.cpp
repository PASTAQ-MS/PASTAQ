#include <algorithm>

#include "metamatch/metamatch.hpp"

void print_metamatch_peaks(const std::vector<MetaMatch::Peak>& peaks) {
    int k = 0;
    for (const auto& peak : peaks) {
        std::cout << "k: " << k;
        std::cout << " mz: " << peak.mz;
        std::cout << " rt: " << peak.rt;
        std::cout << " height: " << peak.height;
        std::cout << " file_id: " << peak.file_id;
        std::cout << " cluster_id: " << peak.cluster_id;
        std::cout << std::endl;
        ++k;
    }
}

void MetaMatch::find_candidates(std::vector<MetaMatch::Peak>& peaks,
                                double radius_mz, double radius_rt) {
    // DEBUG
    // std::cout << "BEFORE:" << std::endl;
    // print_metamatch_peaks(peaks);
    auto sort_peaks = [](auto p1, auto p2) -> bool {
        return (p1.mz < p2.mz) || ((p1.mz == p2.mz) && (p1.rt < p2.rt));
    };
    std::stable_sort(peaks.begin(), peaks.end(), sort_peaks);

    int cluster_id = 0;
    for (size_t i = 0; i < peaks.size(); ++i) {
        auto& peak_a = peaks[i];

        if (peak_a.cluster_id != -1) {
            continue;
        }
        peak_a.cluster_id = cluster_id;

        // Calculate initial centroid stats.
        double x_sum = peak_a.mz * peak_a.height;
        double y_sum = peak_a.rt * peak_a.height;
        double height_sum = peak_a.height;
        double cluster_mz = x_sum / height_sum;
        double cluster_rt = y_sum / height_sum;
        peak_a.cluster_mz = cluster_mz;
        peak_a.cluster_rt = cluster_rt;

        // Mark cluster candidates.
        for (size_t j = (i + 1); j < peaks.size(); ++j) {
            auto& peak_b = peaks[j];
            if (peak_b.mz > cluster_mz + radius_mz &&
                peak_b.rt > cluster_rt + radius_rt) {
                break;
            }
            if (peak_b.cluster_id == -1 && peak_b.file_id != peak_a.file_id &&
                (peak_b.mz < cluster_mz + radius_mz &&
                 peak_b.rt < cluster_rt + radius_rt)) {
                // Add the peak to this cluster.
                peak_b.cluster_id = cluster_id;
                // Update the cluster centroid.
                x_sum += peak_b.mz * peak_b.height;
                y_sum += peak_b.rt * peak_b.height;
                height_sum += peak_b.height;
                cluster_mz = x_sum / height_sum;
                cluster_rt = y_sum / height_sum;
                peak_b.cluster_mz = cluster_mz;
                peak_b.cluster_rt = cluster_rt;
                // Cull far peaks.
                for (size_t k = (i + 1); k <= j; ++k) {
                    auto& peak_c = peaks[k];
                    if (peak_c.cluster_id == cluster_id &&
                        (peak_c.mz > cluster_mz + radius_mz ||
                         peak_c.rt > cluster_rt + radius_rt)) {
                        x_sum -= peak_c.mz * peak_c.height;
                        y_sum -= peak_c.rt * peak_c.height;
                        height_sum -= peak_c.height;
                        cluster_mz = x_sum / height_sum;
                        cluster_rt = y_sum / height_sum;
                        peak_c.cluster_id = -1;
                        peak_c.cluster_mz = peak_c.mz;
                        peak_c.cluster_rt = peak_c.rt;
                    }
                }
            }
        }
        // TODO: Cull multiple peaks per file.
        ++cluster_id;
    }
    // DEBUG
    // std::cout << "AFTER:" << std::endl;
    print_metamatch_peaks(peaks);
    return;
}
