#include <algorithm>
#include <map>
#include <sstream>

#include "metamatch/metamatch.hpp"

std::vector<MetaMatch::FeatureCluster> MetaMatch::find_feature_clusters(
    std::vector<uint64_t>& group_ids,
    std::vector<std::vector<FeatureDetection::Feature>>& features,
    double keep_perc, double intensity_threshold, double n_sig_mz,
    double n_sig_rt) {
    size_t n_files = features.size();

    // We need two sets of indexes, one sorted in descending order of intensity
    // to prioritize the selection of a reference feature to match to, and a set
    // of indexes per file to sort by ascending m/z order. This is done to
    // speedup searching of features using binary search.
    struct Index {
        uint64_t file_index;
        uint64_t group_id;
        uint64_t feature_index;
        double mz;
        double rt;
        double mz_sigma;
        double rt_sigma;
        double intensity;
        int8_t charge_state;
    };

    // Find the total number of features and prepare index vectors.
    auto available_features = std::vector<std::vector<bool>>(n_files);
    auto feature_lists = std::vector<std::vector<Index>>(n_files);
    std::vector<Index> all_features;
    for (size_t i = 0; i < n_files; ++i) {
        size_t n_features = features[i].size();
        available_features[i] = std::vector<bool>(n_features, true);
        feature_lists[i] = std::vector<Index>(n_features);
        for (size_t j = 0; j < n_features; ++j) {
            auto& feature = features[i][j];
            feature_lists[i][j] = {
                i,                                              // file_index
                group_ids[i],                                   // group_id
                j,                                              // feature_index
                feature.monoisotopic_mz,                        // mz
                feature.average_rt + feature.average_rt_delta,  // rt
                feature.average_mz_sigma,                       // mz_sigma
                feature.average_rt_sigma,                       // rt_sigma
                feature.max_height,                             // intensity
                feature.charge_state,                           // charge_state
            };
        }
        // Copy feature_lists[i] to the end of all_features.
        all_features.insert(all_features.end(), feature_lists[i].begin(),
                            feature_lists[i].end());
    }

    // Sort all_features by intensity and feature_lists by mz.
    std::sort(all_features.begin(), all_features.end(),
              [](auto a, auto b) -> bool { return a.intensity > b.intensity; });
    for (auto& feature_list : feature_lists) {
        std::sort(feature_list.begin(), feature_list.end(),
                  [](auto a, auto b) -> bool { return a.mz < b.mz; });
    }

    // Create map with the number of files on each group.
    std::map<uint64_t, uint64_t> groups_map;
    for (const auto& group : group_ids) {
        ++groups_map[group];
    }

    // Start the matching.
    std::vector<MetaMatch::FeatureCluster> clusters;
    size_t cluster_counter = 0;
    for (size_t i = 0; i < all_features.size(); ++i) {
        auto& ref_feature = all_features[i];
        auto ref_file_index = ref_feature.file_index;
        auto ref_feature_index = ref_feature.feature_index;

        // Check availability.
        if (!available_features[ref_file_index][ref_feature_index]) {
            continue;
        }

        // Calculate the boundary region for this feature.
        double ref_min_mz = ref_feature.mz - n_sig_mz * ref_feature.mz_sigma;
        double ref_max_mz = ref_feature.mz + n_sig_mz * ref_feature.mz_sigma;
        double ref_min_rt = ref_feature.rt - n_sig_rt * ref_feature.rt_sigma;
        double ref_max_rt = ref_feature.rt + n_sig_rt * ref_feature.rt_sigma;

        // To search the ROI we use a combination of binary search and linear
        // search. We want to minimize the time we spend on the linear search
        // portion, and for that reason the binary search focuses on m/z, as it
        // is less likely to have points with an exact mass at multiple
        // retention times than the opposite.
        std::vector<MetaMatch::FeatureId> features_in_cluster;
        std::map<uint64_t, uint64_t> cluster_groups_map;
        for (size_t j = 0; j < n_files; ++j) {
            const auto& feature_list = feature_lists[j];
            size_t n_features = feature_list.size();
            size_t left = 0;
            size_t right = n_features;
            while (left < right) {
                size_t mid = left + ((right - left) / 2);
                if (feature_list[mid].mz < ref_min_mz) {
                    left = mid + 1;
                } else {
                    right = mid;
                }
            }
            size_t min_k = right;
            if (right >= n_features || feature_list[min_k].mz > ref_max_mz) {
                continue;
            }

            // Keep track of the best candidate for this file.
            double best_intensity = 0;
            size_t best_index = 0;
            for (size_t k = min_k; k < n_features; ++k) {
                if (feature_list[k].mz > ref_max_mz) {
                    break;
                }
                auto feature = feature_list[k];

                // We are using point-in-rectangle check instead of intersection
                // of boundaries to determine if two features are in range.
                if (feature.charge_state != ref_feature.charge_state ||
                    feature.rt < ref_min_rt || feature.rt > ref_max_rt ||
                    !available_features[j][feature.feature_index]) {
                    continue;
                }

                if (feature.intensity > best_intensity) {
                    best_intensity = feature.intensity;
                    best_index = feature.feature_index;
                }
            }
            if (best_intensity > intensity_threshold) {
                // NOTE: Currently storing the feature index instead of the
                // feature ids for performance. If we are to keep this cluster,
                // they will be swapped.
                features_in_cluster.push_back({j, best_index});
                ++cluster_groups_map[group_ids[j]];
            }
        }

        // Check if the given feature_ids selected for this cluster meet the
        // filter criteria of a minimum number of samples for any given group.
        // For example, if we have three groups of 10 samples, with a
        // `keep_perc' of 0.7, we meet the nan percentage if in any of the three
        // groups we have at least 7 features being matched.
        bool nan_criteria_met = false;
        for (const auto& cluster_group : cluster_groups_map) {
            auto group_id = cluster_group.first;
            auto group_number = cluster_group.second;
            uint64_t required_number = groups_map[group_id] * keep_perc;
            if (group_number >= required_number) {
                nan_criteria_met = true;
                break;
            }
        }
        if (!nan_criteria_met) {
            continue;
        }

        // Build cluster object.
        MetaMatch::FeatureCluster cluster = {};
        cluster.total_heights = std::vector<double>(n_files);
        cluster.monoisotopic_heights = std::vector<double>(n_files);
        cluster.max_heights = std::vector<double>(n_files);
        cluster.total_volumes = std::vector<double>(n_files);
        cluster.monoisotopic_volumes = std::vector<double>(n_files);
        cluster.max_volumes = std::vector<double>(n_files);
        cluster.id = cluster_counter++;
        cluster.charge_state = ref_feature.charge_state;
        for (auto& feature_id : features_in_cluster) {
            size_t file_id = feature_id.file_id;
            size_t feature_index = feature_id.feature_id;
            auto& feature = features[file_id][feature_index];

            // Mark clustered features as not available.
            available_features[file_id][feature_index] = false;

            // Replace feature_index with its corresponding id.
            feature_id.feature_id = feature.id;
            cluster.feature_ids.push_back(feature_id);

            // Store some statistics about the cluster.
            cluster.mz += feature.monoisotopic_mz;
            cluster.rt += feature.average_rt + feature.average_rt_delta;
            cluster.avg_total_height += feature.total_height;
            cluster.avg_monoisotopic_height += feature.monoisotopic_height;
            cluster.avg_max_height += feature.max_height;
            cluster.avg_total_volume += feature.total_volume;
            cluster.avg_monoisotopic_volume += feature.monoisotopic_volume;
            cluster.avg_max_volume += feature.max_volume;
            cluster.total_heights[file_id] = feature.total_height;
            cluster.monoisotopic_heights[file_id] = feature.monoisotopic_height;
            cluster.max_heights[file_id] = feature.max_height;
            cluster.total_volumes[file_id] = feature.total_volume;
            cluster.monoisotopic_volumes[file_id] = feature.monoisotopic_volume;
            cluster.max_volumes[file_id] = feature.max_volume;
        }
        cluster.mz /= features_in_cluster.size();
        cluster.rt /= features_in_cluster.size();
        cluster.avg_total_height /= features_in_cluster.size();
        cluster.avg_monoisotopic_height /= features_in_cluster.size();
        cluster.avg_max_height /= features_in_cluster.size();
        cluster.avg_total_volume /= features_in_cluster.size();
        cluster.avg_monoisotopic_volume /= features_in_cluster.size();
        cluster.avg_max_volume /= features_in_cluster.size();
        clusters.push_back(cluster);
    }

    return clusters;
}

std::vector<MetaMatch::PeakCluster> MetaMatch::find_peak_clusters(
    std::vector<uint64_t>& group_ids,
    std::vector<std::vector<Centroid::Peak>>& peaks,
    double keep_perc, double intensity_threshold, double n_sig_mz,
    double n_sig_rt) {
    size_t n_files = peaks.size();

    // We need two sets of indexes, one sorted in descending order of intensity
    // to prioritize the selection of a reference peak to match to, and a set
    // of indexes per file to sort by ascending m/z order. This is done to
    // speedup searching of peaks using binary search.
    struct Index {
        uint64_t file_index;
        uint64_t group_id;
        uint64_t peak_index;
        double mz;
        double rt;
        double mz_sigma;
        double rt_sigma;
        double intensity;
    };

    // Find the total number of peaks and prepare index vectors.
    auto available_peaks = std::vector<std::vector<bool>>(n_files);
    auto peak_lists = std::vector<std::vector<Index>>(n_files);
    std::vector<Index> all_peaks;
    for (size_t i = 0; i < n_files; ++i) {
        size_t n_peaks = peaks[i].size();
        available_peaks[i] = std::vector<bool>(n_peaks, true);
        peak_lists[i] = std::vector<Index>(n_peaks);
        for (size_t j = 0; j < n_peaks; ++j) {
            auto& peak = peaks[i][j];
            peak_lists[i][j] = {
                i,                                       // file_index
                group_ids[i],                            // group_id
                j,                                       // peak_index
                peak.fitted_mz,                          // mz
                peak.fitted_rt + peak.rt_delta,          // rt
                peak.fitted_sigma_mz,                    // mz_sigma
                peak.fitted_sigma_rt,                    // rt_sigma
                peak.fitted_height,                      // intensity
            };
        }
        // Copy peak_lists[i] to the end of all_peaks.
        all_peaks.insert(all_peaks.end(), peak_lists[i].begin(),
                            peak_lists[i].end());
    }

    // Sort all_peaks by intensity and peak_lists by mz.
    std::sort(all_peaks.begin(), all_peaks.end(),
              [](auto a, auto b) -> bool { return a.intensity > b.intensity; });
    for (auto& peak_list : peak_lists) {
        std::sort(peak_list.begin(), peak_list.end(),
                  [](auto a, auto b) -> bool { return a.mz < b.mz; });
    }

    // Create map with the number of files on each group.
    std::map<uint64_t, uint64_t> groups_map;
    for (const auto& group : group_ids) {
        ++groups_map[group];
    }

    // Start the matching.
    std::vector<MetaMatch::PeakCluster> clusters;
    size_t cluster_counter = 0;
    for (size_t i = 0; i < all_peaks.size(); ++i) {
        auto& ref_peak = all_peaks[i];
        auto ref_file_index = ref_peak.file_index;
        auto ref_peak_index = ref_peak.peak_index;

        // Check availability.
        if (!available_peaks[ref_file_index][ref_peak_index]) {
            continue;
        }

        // Calculate the boundary region for this peak.
        double ref_min_mz = ref_peak.mz - n_sig_mz * ref_peak.mz_sigma;
        double ref_max_mz = ref_peak.mz + n_sig_mz * ref_peak.mz_sigma;
        double ref_min_rt = ref_peak.rt - n_sig_rt * ref_peak.rt_sigma;
        double ref_max_rt = ref_peak.rt + n_sig_rt * ref_peak.rt_sigma;

        // To search the ROI we use a combination of binary search and linear
        // search. We want to minimize the time we spend on the linear search
        // portion, and for that reason the binary search focuses on m/z, as it
        // is less likely to have points with an exact mass at multiple
        // retention times than the opposite.
        std::vector<MetaMatch::PeakId> peaks_in_cluster;
        std::map<uint64_t, uint64_t> cluster_groups_map;
        for (size_t j = 0; j < n_files; ++j) {
            const auto& peak_list = peak_lists[j];
            size_t n_peaks = peak_list.size();
            size_t left = 0;
            size_t right = n_peaks;
            while (left < right) {
                size_t mid = left + ((right - left) / 2);
                if (peak_list[mid].mz < ref_min_mz) {
                    left = mid + 1;
                } else {
                    right = mid;
                }
            }
            size_t min_k = right;
            if (right >= n_peaks || peak_list[min_k].mz > ref_max_mz) {
                continue;
            }

            // Keep track of the best candidate for this file.
            double best_intensity = 0;
            size_t best_index = 0;
            for (size_t k = min_k; k < n_peaks; ++k) {
                if (peak_list[k].mz > ref_max_mz) {
                    break;
                }
                auto peak = peak_list[k];

                // We are using point-in-rectangle check instead of intersection
                // of boundaries to determine if two peaks are in range.
                if (peak.rt < ref_min_rt || peak.rt > ref_max_rt ||
                    !available_peaks[j][peak.peak_index]) {
                    continue;
                }

                if (peak.intensity > best_intensity) {
                    best_intensity = peak.intensity;
                    best_index = peak.peak_index;
                }
            }
            if (best_intensity > intensity_threshold) {
                // NOTE: Currently storing the peak index instead of the
                // peak ids for performance. If we are to keep this cluster,
                // they will be swapped.
                peaks_in_cluster.push_back({j, best_index});
                ++cluster_groups_map[group_ids[j]];
            }
        }

        // Check if the given peak_ids selected for this cluster meet the
        // filter criteria of a minimum number of samples for any given group.
        // For example, if we have three groups of 10 samples, with a
        // `keep_perc' of 0.7, we meet the nan percentage if in any of the three
        // groups we have at least 7 peaks being matched.
        bool nan_criteria_met = false;
        for (const auto& cluster_group : cluster_groups_map) {
            auto group_id = cluster_group.first;
            auto group_number = cluster_group.second;
            uint64_t required_number = groups_map[group_id] * keep_perc;
            if (group_number >= required_number) {
                nan_criteria_met = true;
                break;
            }
        }
        if (!nan_criteria_met) {
            continue;
        }

        // Build cluster object.
        MetaMatch::PeakCluster cluster = {};
        cluster.heights = std::vector<double>(n_files);
        cluster.volumes = std::vector<double>(n_files);
        cluster.id = cluster_counter++;
        for (auto& peak_id : peaks_in_cluster) {
            size_t file_id = peak_id.file_id;
            size_t peak_index = peak_id.peak_id;
            auto& peak = peaks[file_id][peak_index];

            // Mark clustered peaks as not available.
            available_peaks[file_id][peak_index] = false;

            // Replace peak_index with its corresponding id.
            peak_id.peak_id = peak.id;
            cluster.peak_ids.push_back(peak_id);

            // Store some statistics about the cluster.
            cluster.mz += peak.fitted_mz;
            cluster.rt += peak.fitted_rt + peak.rt_delta;
            cluster.avg_height += peak.fitted_height;
            cluster.avg_volume += peak.fitted_volume;
            cluster.heights[file_id] = peak.fitted_height;
            cluster.volumes[file_id] = peak.fitted_volume;
        }
        cluster.mz /= peaks_in_cluster.size();
        cluster.rt /= peaks_in_cluster.size();
        cluster.avg_height /= peaks_in_cluster.size();
        cluster.avg_volume /= peaks_in_cluster.size();
        clusters.push_back(cluster);
    }

    return clusters;
}
