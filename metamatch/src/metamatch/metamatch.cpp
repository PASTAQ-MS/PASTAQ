#include <algorithm>
#include <sstream>

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
void update_centroid(MetaMatch::Peak& peak, int cluster_id, double& x_sum,
                     double& y_sum, double& height_sum, double& cluster_mz,
                     double& cluster_rt) {
    // Add the peak to this cluster.
    peak.cluster_id = cluster_id;
    // Update the cluster centroid.
    x_sum += peak.mz * peak.height;
    y_sum += peak.rt * peak.height;
    height_sum += peak.height;
    cluster_mz = x_sum / height_sum;
    cluster_rt = y_sum / height_sum;
    peak.cluster_mz = cluster_mz;
    peak.cluster_rt = cluster_rt;
}

void cull_far_peaks(std::vector<MetaMatch::Peak>& peaks, int cluster_id,
                    size_t i, size_t j, double radius_mz, double radius_rt,
                    double& x_sum, double& y_sum, double& height_sum,
                    double& cluster_mz, double& cluster_rt) {
    // Cull far peaks.
    for (size_t k = (i + 1); k <= j; ++k) {
        auto& peak = peaks[k];
        if (peak.cluster_id == cluster_id &&
            (peak.mz > cluster_mz + radius_mz ||
             peak.rt > cluster_rt + radius_rt)) {
            x_sum -= peak.mz * peak.height;
            y_sum -= peak.rt * peak.height;
            height_sum -= peak.height;
            cluster_mz = x_sum / height_sum;
            cluster_rt = y_sum / height_sum;
            peak.cluster_id = -1;
            peak.cluster_mz = peak.mz;
            peak.cluster_rt = peak.rt;
        }
    }
}

bool check_fraction(std::vector<MetaMatch::Peak>& peaks, int cluster_id,
                    size_t i, double cluster_mz, double cluster_rt,
                    const MetaMatch::Parameters& parameters) {
    // Check if this cluster contains the
    // necessary number of peaks per class. If the
    // fraction is not big enough, free the peaks
    // for future clustering.
    std::vector<size_t> file_ids;
    std::vector<size_t> class_map(parameters.n_classes);
    for (size_t j = i; j < peaks.size(); ++j) {
        auto& peak = peaks[j];
        if (peak.mz > cluster_mz + parameters.radius_mz) {
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
        if ((double)n_hits / (double)parameters.n_files > parameters.fraction) {
            fraction_achieved = true;
            break;
        }
    }
    return fraction_achieved;
}

void MetaMatch::find_candidates(std::vector<MetaMatch::Peak>& peaks,
                                const MetaMatch::Parameters& parameters) {
    // DEBUG
    // std::cout << "BEFORE:" << std::endl;
    // print_metamatch_peaks(peaks);
    std::cout << "sorting peaks..." << std::endl;
    auto sort_peaks = [](auto p1, auto p2) -> bool {
        return (p1.mz < p2.mz) || ((p1.mz == p2.mz) && (p1.rt < p2.rt)) ||
               ((p1.rt == p2.rt) && (p1.file_id < p2.file_id));
    };
    std::stable_sort(peaks.begin(), peaks.end(), sort_peaks);
    // DEBUG
    // print_metamatch_peaks(peaks);

    std::cout << "clustering..." << std::endl;
    int cluster_id = 0;
    double avg_peaks_per_iter = 0;
    for (size_t i = 0; i < peaks.size(); ++i) {
        auto& peak_a = peaks[i];
        // DEBUG
        if (i % 10000 == 0) {
            std::cout << "progress: peak " << i << " out of " << peaks.size()
                      << std::endl;
        }

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
        size_t visited = 1;
        for (size_t j = (i + 1); j < peaks.size(); ++j, ++visited) {
            auto& peak_b = peaks[j];
            // Since we know that the peaks are sorted monotonically in mz and
            // then rt, in order to calculate the maximum potential j we only
            // need to find the point where the peak.mz is above the cluster
            // radius.
            if (peak_b.mz > cluster_mz + parameters.radius_mz) {
                avg_peaks_per_iter += visited;
                break;
            }
            if (peak_b.cluster_id == -1 && peak_b.file_id != peak_a.file_id &&
                (peak_b.mz < cluster_mz + parameters.radius_mz &&
                 peak_b.rt < cluster_rt + parameters.radius_rt)) {
                // FIXME: update_centroid(peak_b, cluster_id, x_sum, y_sum,
                // height_sum, cluster_mz, cluster_rt); Add the peak to this
                // cluster.
                peak_b.cluster_id = cluster_id;
                // Update the cluster centroid.
                x_sum += peak_b.mz * peak_b.height;
                y_sum += peak_b.rt * peak_b.height;
                height_sum += peak_b.height;
                cluster_mz = x_sum / height_sum;
                cluster_rt = y_sum / height_sum;
                peak_b.cluster_mz = cluster_mz;
                peak_b.cluster_rt = cluster_rt;

                // FIXME: cull_far_peaks(peaks, cluster_id, i, j,
                // parameters.radius_mz,
                // parameters.radius_rt, x_sum, y_sum, height_sum,
                // cluster_mz, cluster_rt);
                // Cull far peaks.
                for (size_t k = (i + 1); k <= j; ++k) {
                    auto& peak = peaks[k];
                    if (peak.cluster_id == cluster_id &&
                        (peak.mz > cluster_mz + parameters.radius_mz ||
                         peak.rt > cluster_rt + parameters.radius_rt)) {
                        x_sum -= peak.mz * peak.height;
                        y_sum -= peak.rt * peak.height;
                        height_sum -= peak.height;
                        cluster_mz = x_sum / height_sum;
                        cluster_rt = y_sum / height_sum;
                        peak.cluster_id = -1;
                        peak.cluster_mz = peak.mz;
                        peak.cluster_rt = peak.rt;
                    }
                }
            }
        }
        // TODO: Cull multiple peaks per file.
        // ...

        if (!check_fraction(peaks, cluster_id, i, cluster_mz, cluster_rt,
                            parameters)) {
            for (size_t j = i; j < peaks.size(); ++j) {
                auto& peak = peaks[j];
                if (peak.mz > cluster_mz + parameters.radius_mz) {
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
    std::cout << "SUM PEAKS VISITED: " << avg_peaks_per_iter << std::endl;
    std::cout << "AVG PEAKS VISITED PER ITER: "
              << avg_peaks_per_iter / (double)peaks.size() << std::endl;

    // DEBUG
    // print_metamatch_peaks(peaks);
    return;
}

std::vector<MetaMatch::Peak> MetaMatch::extract_orphans(
    std::vector<MetaMatch::Peak>& peaks) {
    std::vector<MetaMatch::Peak> orphans;
    // TODO(alex): Where does the sorting need to
    // happen?
    auto sort_by_cluster_id = [](auto p1, auto p2) -> bool {
        return p1.cluster_id < p2.cluster_id;
    };
    std::stable_sort(peaks.begin(), peaks.end(), sort_by_cluster_id);

    for (size_t i = 0; i < peaks.size(); ++i) {
        if (peaks[i].cluster_id != -1) {
            orphans.insert(orphans.end(), peaks.begin(), peaks.begin() + i);
            peaks.erase(peaks.begin(), peaks.begin() + i);
            break;
        }
    }
    // DEBUG
    // std::cout << "ORPHANS:" << std::endl;
    // print_metamatch_peaks(orphans);
    // std::cout << "NOT ORPHANS:" << std::endl;
    // print_metamatch_peaks(peaks);
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
    std::stable_sort(peaks.begin(), peaks.end(), sort_by_cluster_id);
    if (peaks.empty()) {
        return clusters;
    }
    // DEBUG
    // std::cout << "SORTED METAPEAKS:" << std::endl;
    // print_metamatch_peaks(peaks);

    int cluster_id = peaks[0].cluster_id;
    size_t k = 0;
    for (size_t i = 1; i < peaks.size(); ++i) {
        if (cluster_id != peaks[i].cluster_id) {
            MetaMatch::Cluster cluster = {};
            cluster.id = cluster_id;
            cluster.mz = peaks[i - 1].cluster_mz;
            cluster.rt = peaks[i - 1].cluster_rt;
            cluster.file_heights = std::vector<double>(n_files);
            for (size_t file_index = 0; file_index < n_files; ++file_index) {
                if (file_index == peaks[k].file_id) {
                    cluster.file_heights[file_index] = peaks[k].height;
                    ++k;
                } else {
                    cluster.file_heights[file_index] = 0;
                }
            }
            clusters.push_back(cluster);
            cluster_id = peaks[i].cluster_id;
            k = i;
        }
    }

    // DEBUG
    // Print clusters...
    // for (const auto& cluster : clusters) {
    // std::cout << "id: " << cluster.id;
    // std::cout << " mz: " << cluster.mz;
    // std::cout << " rt: " << cluster.rt;
    // std::cout << " elements: [";
    // for (const auto& height : cluster.file_heights)
    // { std::cout << height << ",";
    //}
    // std::cout << "]" << std::endl;
    //}
    return clusters;
}

// TODO(alex): Move to MetaMatch::Files::Csv::write_clusters
bool MetaMatch::write_clusters(std::ostream& stream,
                               const std::vector<MetaMatch::Cluster>& clusters,
                               const MetaMatch::Parameters& parameters) {
    char cell_delimiter = ' ';
    char line_delimiter = '\n';

    // Prepare the CSV header.
    std::vector<std::string> header_columns = {
        "metapeak",
        "mz",
        "rt",
    };
    for (size_t i = 0; i < parameters.n_files; ++i) {
        header_columns.push_back("file_h" + std::to_string(i));
    }
    // Write the CSV header.
    for (size_t i = 0; i < header_columns.size(); ++i) {
        stream << header_columns[i];
        if (i == header_columns.size() - 1) {
            stream << line_delimiter;
        } else {
            stream << cell_delimiter;
        }
    }
    stream.precision(8);
    for (const auto& cluster : clusters) {
        // metapeak
        stream << cluster.id << cell_delimiter;
        // mz
        stream << cluster.mz << cell_delimiter;
        // rt
        stream << cluster.rt << cell_delimiter;
        // file heights
        for (size_t i = 0; i < parameters.n_files; ++i) {
            stream << cluster.file_heights[i];
            if (i == parameters.n_files - 1) {
                stream << line_delimiter;
            } else {
                stream << cell_delimiter;
            }
        }
    }
    return stream.good();
}

std::vector<std::pair<std::string, size_t>> MetaMatch::read_file_list(
    std::istream& stream) {
    char cell_delimiter = ' ';
    char line_delimiter = '\n';
    std::vector<std::pair<std::string, size_t>> ret;
    std::string line;
    // TODO(alex): Add validation (return bool).
    while (std::getline(stream, line, line_delimiter)) {
        std::string token;
        std::stringstream token_stream(line);
        std::pair<std::string, size_t> file = {};
        // TODO: validation...
        std::getline(token_stream, token, cell_delimiter);
        std::istringstream(token) >> file.first;
        // TODO: validation...
        std::getline(token_stream, token, cell_delimiter);
        std::istringstream(token) >> file.second;
        ret.push_back(file);
    }
    return ret;
}
