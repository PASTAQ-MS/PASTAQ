#include <algorithm>

#include "metamatch/metamatch.hpp"

void MetaMatch::find_candidates(
    const std::vector<std::vector<Centroid::Peak>>& peak_files) {
    for (size_t file_index = 0; file_index < peak_files.size(); ++file_index) {
        const auto& file_a = peak_files[file_index];
        // DEBUG
        std::cout << "file: " << file_index << std::endl;
        for (size_t peak_index = 0; peak_index < file_a.size(); ++peak_index) {
            const auto& peak_a = file_a[peak_index];
            // Calculate peak_a region.
            double min_mz_a = peak_a.mz - 3 * peak_a.sigma_mz;
            double max_mz_a = peak_a.mz + 3 * peak_a.sigma_mz;
            double min_rt_a = peak_a.rt - 3 * peak_a.sigma_rt;
            double max_rt_a = peak_a.rt + 3 * peak_a.sigma_rt;

            // DEBUG
            std::cout << "peak [" << peak_index << "]:";
            std::cout << " mz: " << peak_a.mz;
            std::cout << " rt: " << peak_a.rt;
            std::cout << std::endl;
            // Find candidate peaks.
            // FIXME: Copying entire peaks multiple times. Memory and
            // performance limitations incoming! We really want to store
            // pointers or indices instead but implementing this first like this
            // for simplicity in the exploration.
            std::vector<Centroid::Peak> candidates;
            for (size_t i = 0; i < peak_files.size(); ++i) {
                if (i == file_index) {
                    continue;
                }
                // Find peaks in +/- sigma region.
                const auto& file_b = peak_files[i];
                for (size_t j = 0; j < file_b.size(); ++j) {
                    const auto& peak_b = file_b[j];
                    // Calculate peak_b region.
                    double min_mz_b = peak_b.mz - 3 * peak_b.sigma_mz;
                    double max_mz_b = peak_b.mz + 3 * peak_b.sigma_mz;
                    double min_rt_b = peak_b.rt - 3 * peak_b.sigma_rt;
                    double max_rt_b = peak_b.rt + 3 * peak_b.sigma_rt;
                    if (max_rt_a < min_rt_b || max_rt_b < min_rt_a ||
                        max_mz_a < min_mz_b || max_mz_b < min_mz_a) {
                        continue;
                    }
                    candidates.push_back(peak_b);
                }
            }
            // DEBUG
            int k = 0;
            for (const auto& peak_b : candidates) {
                std::cout << "candidate [" << k << "]:";
                std::cout << " mz: " << peak_b.mz;
                std::cout << " rt: " << peak_b.rt;
                std::cout << std::endl;
                ++k;
            }
        }
    }
    return;
}

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
    // print_metamatch_peaks(peaks);
    return;
}
