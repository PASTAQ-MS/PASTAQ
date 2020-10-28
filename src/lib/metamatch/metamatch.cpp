#include <algorithm>
#include <map>
#include <sstream>

#include "metamatch/metamatch.hpp"

void calculate_cluster_pos(double& cluster_mz, double& cluster_rt,
                           double& cluster_sigma_mz, double& cluster_sigma_rt,
                           std::vector<MetaMatch::Peak>& peaks,
                           std::vector<size_t>& metapeak_indexes) {
    double avg_mz = 0;
    double avg_rt = 0;
    double avg_sigma_mz = 0;
    double avg_sigma_rt = 0;
    double height_sum = 0;
    for (const auto& index : metapeak_indexes) {
        double mz = peaks[index].fitted_mz;
        double rt = peaks[index].fitted_rt + peaks[index].rt_delta;
        double sigma_mz = peaks[index].fitted_sigma_mz;
        double sigma_rt = peaks[index].fitted_sigma_rt;
        double height = peaks[index].fitted_height;
        avg_mz += mz * height;
        avg_rt += rt * height;
        avg_sigma_mz += sigma_mz * height;
        avg_sigma_rt += sigma_rt * height;
        height_sum += height;

        // NOTE(alex): This is a weighted average, it might cause problems if
        // the overall intensity between the files are very different (The
        // centroid will be biased towards the file with the greater average
        // intensity). If instead of a weighted average we want a density
        // average we could use the following:
        //
        //     avg_mz += peaks[index].fitted_mz;
        //     avg_rt += peaks[index].fitted_rt;
        //     height_sum += 1;
        //
        // This might have a problem where the noise could have a greater impact
        // in the position of the centroid.
    }
    cluster_mz = avg_mz / height_sum;
    cluster_rt = avg_rt / height_sum;
    cluster_sigma_mz = avg_sigma_mz / height_sum;
    cluster_sigma_rt = avg_sigma_rt / height_sum;
}

void MetaMatch::find_clusters(std::vector<MetaMatch::Peak>& peaks,
        const std::vector<ClassMap>& class_maps, double n_sig_mz, double n_sig_rt) {
    auto sort_peaks = [](auto p1, auto p2) -> bool {
        double p1_mz = p1.fitted_mz;
        double p2_mz = p2.fitted_mz;
        double p1_rt = p1.fitted_rt + p1.rt_delta;
        double p2_rt = p2.fitted_rt + p2.rt_delta;
        return (p1_mz < p2_mz) || ((p1_mz == p2_mz) && (p1_rt < p2_rt)) ||
               ((p1_rt == p2_rt) && (p1.file_id < p2.file_id));
    };
    std::sort(peaks.begin(), peaks.end(), sort_peaks);

    int cluster_id = 0;
    for (size_t i = 0; i < peaks.size(); ++i) {
        if (peaks[i].cluster_id != -1) {
            continue;
        }

        // Calculate initial centroid stats.
        double cluster_mz = 0;
        double cluster_rt = 0;
        double cluster_sigma_mz = 0;
        double cluster_sigma_rt = 0;
        std::vector<size_t> metapeak_indexes = {i};
        peaks[i].cluster_id = cluster_id;
        calculate_cluster_pos(cluster_mz, cluster_rt, cluster_sigma_mz,
                              cluster_sigma_rt, peaks, metapeak_indexes);

        for (size_t j = (i + 1); j < peaks.size(); ++j) {
            // Since we know that the peaks are sorted monotonically in mz and
            // then rt, in order to calculate the maximum potential j we only
            // need to find the point where the peak.fitted_mz is above the
            // cluster radius.
            if (peaks[j].fitted_mz > cluster_mz + cluster_sigma_mz) {
                break;
            }
            double peak_mz = peaks[j].fitted_mz;
            double peak_rt = peaks[j].fitted_rt + peaks[j].rt_delta;
            if (peaks[j].cluster_id == -1 &&
                (peak_mz < cluster_mz + cluster_sigma_mz &&
                 peak_rt < cluster_rt + cluster_sigma_rt)) {
                // If the cluster already contains a peak from the same file as
                // peaks[j], check if height of said peak is greater than
                // peaks[j].fitted_height, if it is, swap the index,
                // otherwise, continue.
                bool file_found = false;
                for (auto& index : metapeak_indexes) {
                    if (peaks[index].file_id == peaks[j].file_id) {
                        file_found = true;
                        if (peaks[index].file_id == peaks[j].file_id &&
                            peaks[index].fitted_height <
                                peaks[j].fitted_height) {
                            // Update cluster peaks.
                            peaks[index].cluster_id = -1;
                            index = j;
                            peaks[index].cluster_id = cluster_id;
                            break;
                        }
                    }
                }
                if (!file_found) {
                    peaks[j].cluster_id = cluster_id;
                    metapeak_indexes.push_back(j);
                }
                calculate_cluster_pos(cluster_mz, cluster_rt, cluster_sigma_mz,
                                      cluster_sigma_rt, peaks,
                                      metapeak_indexes);

                // Cull far peaks.
                for (int k = metapeak_indexes.size() - 1; k >= 0; --k) {
                    auto& index = metapeak_indexes[k];
                    double peak_mz = peaks[index].fitted_mz;
                    double peak_rt =
                        peaks[index].fitted_rt + peaks[index].rt_delta;
                    if (peak_mz > cluster_mz + n_sig_mz * cluster_sigma_mz ||
                        peak_mz < cluster_mz - n_sig_mz * cluster_sigma_mz ||
                        peak_rt > cluster_rt + n_sig_rt * cluster_sigma_rt ||
                        peak_rt < cluster_rt - n_sig_rt * cluster_sigma_rt) {
                        peaks[index].cluster_id = -1;
                        metapeak_indexes.erase(metapeak_indexes.begin() + k);
                        calculate_cluster_pos(
                            cluster_mz, cluster_rt, cluster_sigma_mz,
                            cluster_sigma_rt, peaks, metapeak_indexes);
                    }
                }
            }
        }

        // Check if there is enough hits from a class in order to consider the
        // cluster valid.
        std::vector<size_t> class_hits(class_maps.size());
        for (const auto& index : metapeak_indexes) {
            const auto& peak = peaks[index];
            for (size_t k = 0; k < class_hits.size(); ++k) {
                if (peak.class_id == class_maps[k].id) {
                    class_hits[k] += 1;
                    break;
                }
            }
        }
        bool fraction_achieved = false;
        for (size_t k = 0; k < class_hits.size(); ++k) {
            const auto& n_hits = class_hits[k];
            if (n_hits >= class_maps[k].required_hits && n_hits != 0) {
                fraction_achieved = true;
                break;
            }
        }
        if (!fraction_achieved) {
            for (const auto& index : metapeak_indexes) {
                peaks[index].cluster_id = -1;
            }
        } else {
            for (const auto& index : metapeak_indexes) {
                peaks[index].cluster_mz = cluster_mz;
                peaks[index].cluster_rt = cluster_rt;
            }
            ++cluster_id;
        }
    }
}

std::vector<MetaMatch::Peak> MetaMatch::extract_orphans(
    std::vector<MetaMatch::Peak>& peaks) {
    std::vector<MetaMatch::Peak> orphans;
    // TODO(alex): Where does the sorting need to
    // happen?
    auto sort_by_cluster_id = [](auto p1, auto p2) -> bool {
        return p1.cluster_id < p2.cluster_id;
    };
    std::sort(peaks.begin(), peaks.end(), sort_by_cluster_id);

    for (size_t i = 0; i < peaks.size(); ++i) {
        if (peaks[i].cluster_id != -1) {
            orphans.insert(orphans.end(), peaks.begin(), peaks.begin() + i);
            peaks.erase(peaks.begin(), peaks.begin() + i);
            break;
        }
    }
    return orphans;
}

std::vector<MetaMatch::Cluster> MetaMatch::reduce_cluster(
    std::vector<MetaMatch::Peak>& peaks, size_t n_files) {
    std::vector<MetaMatch::Cluster> clusters;
    // TODO(alex): Where does the sorting need to
    // happen?
    auto sort_by_cluster_id = [](auto p1, auto p2) -> bool {
        return p1.cluster_id < p2.cluster_id ||
               (p1.cluster_id == p2.cluster_id && p1.file_id < p2.file_id);
    };
    std::sort(peaks.begin(), peaks.end(), sort_by_cluster_id);
    if (peaks.empty()) {
        return clusters;
    }

    int64_t cluster_id = 0;
    size_t k = 0;
    for (size_t i = 1; i < peaks.size(); ++i) {
        if (cluster_id != peaks[i].cluster_id) {
            MetaMatch::Cluster cluster = {};
            cluster.id = cluster_id;
            cluster.mz = peaks[i - 1].cluster_mz;
            cluster.rt = peaks[i - 1].cluster_rt;
            cluster.file_heights = std::vector<double>(n_files);
            cluster.file_volumes = std::vector<double>(n_files);
            double sum_height = 0.0;
            size_t hits = 0;
            for (size_t j = k; j < i; ++j) {
                auto file_id = peaks[j].file_id;
                cluster.file_heights[file_id] = peaks[j].fitted_height;
                cluster.file_volumes[file_id] = peaks[j].fitted_volume;
                sum_height += peaks[j].fitted_height;
                hits++;
                cluster.peak_ids.push_back({file_id, peaks[j].id});
            }
            cluster.avg_height = sum_height / hits;
            clusters.push_back(cluster);
            cluster_id = peaks[i].cluster_id;
            k = i;
        }
    }

    // Sort clusters by highest average height and change cluster ids
    // accordingly.
    auto sort_by_avg_height = [](auto c1, auto c2) -> bool {
        return c1.avg_height > c2.avg_height;
    };
    std::sort(clusters.begin(), clusters.end(), sort_by_avg_height);
    size_t i = 0;
    for (auto& cluster : clusters) {
        cluster.id = i;
        i++;
    }

    return clusters;
}

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
