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
                                double radius_mz, double radius_rt,
                                size_t n_files, size_t n_classes,
                                double fraction) {
    // DEBUG
    // std::cout << "BEFORE:" << std::endl;
    // print_metamatch_peaks(peaks);
    auto sort_peaks = [](auto p1, auto p2) -> bool {
        return (p1.mz < p2.mz) || ((p1.mz == p2.mz) && (p1.rt < p2.rt)) ||
               ((p1.rt == p2.rt) && (p1.file_id < p2.file_id));
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
        // ...

        // Check if this cluster contains the necessary number of peaks
        // per class. If the fraction is not big enough, free the peaks for
        // future clustering.
        std::vector<size_t> file_ids;
        std::vector<size_t> class_map(n_classes);
        for (size_t j = i; j < peaks.size(); ++j) {
            auto& peak = peaks[j];
            if (peak.mz > cluster_mz + radius_mz &&
                peak.rt > cluster_rt + radius_rt) {
                break;
            }
            if (peak.cluster_id == cluster_id) {
                bool id_present = false;
                for (const auto& id : file_ids) {
                    if (id == peak.file_id) {
                        id_present = true;
                        break;
                    }
                }
                if (!id_present) {
                    file_ids.push_back(peak.file_id);
                    class_map[peak.class_id] += 1;
                }
            }
        }
        bool fraction_achieved = false;
        for (const auto& n_hits : class_map) {
            if ((double)n_hits / (double)n_files > fraction) {
                fraction_achieved = true;
                break;
            }
        }
        if (!fraction_achieved) {
            for (size_t j = i; j < peaks.size(); ++j) {
                auto& peak = peaks[j];
                if (peak.mz > cluster_mz + radius_mz &&
                    peak.rt > cluster_rt + radius_rt) {
                    break;
                }
                if (peak.cluster_id == cluster_id) {
                    peak.cluster_id = -1;
                    peak.cluster_mz = peak.mz;
                    peak.cluster_rt = peak.rt;
                }
            }
        } else {
            ++cluster_id;
        }
    }

    // DEBUG
    print_metamatch_peaks(peaks);
    return;
}

std::vector<MetaMatch::Peak> MetaMatch::extract_orphans(
    std::vector<MetaMatch::Peak>& peaks) {
    std::vector<MetaMatch::Peak> orphans;
    // TODO(alex): Where does the sorting need to happen?
    auto sort_by_cluster_id = [](auto p1, auto p2) -> bool {
        return p1.cluster_id < p2.cluster_id;
    };
    std::stable_sort(peaks.begin(), peaks.end(), sort_by_cluster_id);

    for (size_t i = 0; i < peaks.size(); ++i) {
        if (peaks[i].cluster_id != -1) {
            orphans.insert(orphans.end(),
                           std::make_move_iterator(peaks.begin()),
                           std::make_move_iterator(peaks.begin() + i));
            peaks.erase(peaks.begin(), peaks.begin() + i);
            break;
        }
    }
    // DEBUG
    std::cout << "ORPHANS:" << std::endl;
    print_metamatch_peaks(orphans);
    std::cout << "NOT ORPHANS:" << std::endl;
    print_metamatch_peaks(peaks);
    return orphans;
}
