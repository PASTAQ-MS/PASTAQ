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
        std::cout << " class_id: " << peak.class_id;
        std::cout << " cluster_id: " << peak.cluster_id;
        std::cout << std::endl;
        ++k;
    }
}

void MetaMatch::find_candidates(std::vector<MetaMatch::Peak>& peaks,
                                const MetaMatch::Parameters& parameters) {
    // DEBUG
    // std::cout << "BEFORE:" << std::endl;
    // print_metamatch_peaks(peaks);
    std::cout << "Sorting peaks..." << std::endl;
    auto sort_peaks = [](auto p1, auto p2) -> bool {
        return (p1.mz < p2.mz) || ((p1.mz == p2.mz) && (p1.rt < p2.rt)) ||
               ((p1.rt == p2.rt) && (p1.file_id < p2.file_id));
    };
    std::stable_sort(peaks.begin(), peaks.end(), sort_peaks);
    // DEBUG
    // print_metamatch_peaks(peaks);
    // TODO(alex): Do a first pass through the peaks to calculate the number of
    // files, number of classes, number of files per class.
    std::cout << "Calculating class/file proportion..." << std::endl;
    std::vector<size_t> file_ids;
    struct ClassMeta {
        size_t id;
        std::vector<size_t> file_ids;
    };
    std::vector<ClassMeta> classes;
    for (const auto& peak : peaks) {
        bool class_found = false;
        for (size_t i = 0; i < classes.size(); ++i) {
            if (peak.class_id == classes[i].id) {
                class_found = true;
                bool file_found = false;
                for (size_t j = 0; j < classes[i].file_ids.size(); ++j) {
                    if (peak.file_id == classes[i].file_ids[j]) {
                        file_found = true;
                        break;
                    }
                }
                if (!file_found) {
                    classes[i].file_ids.push_back(peak.file_id);
                }
                break;
            }
        }
        if (!class_found) {
            classes.push_back({peak.class_id, {peak.file_id}});
        }
    }
    for (const auto& cls : classes) {
        std::cout << "class_id: " << cls.id
                  << " n_files: " << cls.file_ids.size() << std::endl;
    }

    std::cout << "Clustering..." << std::endl;
    int cluster_id = 0;
    double avg_peaks_per_iter = 0;
    for (size_t i = 0; i < peaks.size(); ++i) {
        // DEBUG
        // std::cout << "i: " << i << std::endl;
        // print_metamatch_peaks(peaks);
        auto& peak_a = peaks[i];
        // DEBUG
        if (i % 100000 == 0) {
            std::cout << "Progress: peak " << i << " out of " << peaks.size()
                      << std::endl;
        }

        if (peak_a.cluster_id != -1) {
            continue;
        }

        // Calculate initial centroid stats.
        double cluster_mz = 0;
        double cluster_rt = 0;
        std::vector<size_t> metapeak_indexes = {i};
        auto calculate_cluster_pos = [&cluster_mz, &cluster_rt, &peaks,
                                      &metapeak_indexes]() {
            double x_sum = 0;
            double y_sum = 0;
            double height_sum = 0;
            for (const auto& index : metapeak_indexes) {
                x_sum += peaks[index].mz * peaks[index].height;
                y_sum += peaks[index].rt * peaks[index].height;
                height_sum += peaks[index].height;
            }
            cluster_mz = x_sum / height_sum;
            cluster_rt = y_sum / height_sum;
        };
        // auto calculate_cluster_pos = [&cluster_mz, &cluster_rt, &peaks,
        //&metapeak_indexes]() {
        // double x_sum = 0;
        // double y_sum = 0;
        // double height_sum = 0;
        // for (const auto& index : metapeak_indexes) {
        // x_sum += peaks[index].mz;
        // y_sum += peaks[index].rt;
        // height_sum += 1;
        //}
        // cluster_mz = x_sum / height_sum;
        // cluster_rt = y_sum / height_sum;
        //};
        calculate_cluster_pos();

        peak_a.cluster_id = cluster_id;

        // Mark cluster candidates.
        size_t visited = 1;  // DEBUG
        for (size_t j = (i + 1); j < peaks.size(); ++j, ++visited) {
            auto& peak_b = peaks[j];
            // Since we know that the peaks are sorted monotonically in mz and
            // then rt, in order to calculate the maximum potential j we only
            // need to find the point where the peak.mz is above the cluster
            // radius.
            if (peak_b.mz > cluster_mz + parameters.radius_mz) {
                avg_peaks_per_iter += visited;  // DEBUG
                break;
            }
            if (peak_b.cluster_id == -1 && peak_b.file_id != peak_a.file_id &&
                (peak_b.mz < cluster_mz + parameters.radius_mz &&
                 peak_b.rt < cluster_rt + parameters.radius_rt)) {
                // If the cluster already contains a peak from the same file as
                // peak_b, check if height of said peak is greater than
                // peak_b.height, if it is, swap the index, otherwise, continue.
                bool file_found = false;
                for (auto& index : metapeak_indexes) {
                    if (peaks[index].file_id == peak_b.file_id &&
                        peaks[index].height < peak_b.height) {
                        // Update cluster peaks.
                        peaks[index].cluster_id = -1;
                        index = j;
                        peaks[index].cluster_id = cluster_id;
                        calculate_cluster_pos();
                        file_found = true;
                        break;
                    }
                }
                if (!file_found) {
                    peak_b.cluster_id = cluster_id;
                    metapeak_indexes.push_back(j);
                    calculate_cluster_pos();
                }
                // DEBUG
                // std::cout << "i: " << i
                //<< " peaks in cluster: " << metapeak_indexes.size()
                //<< " cluster_mz: " << cluster_mz
                //<< " cluster_rt: " << cluster_rt << std::endl;
                // Cull far peaks.
                for (int k = metapeak_indexes.size() - 1; k >= 0; --k) {
                    auto& index = metapeak_indexes[k];
                    if (peaks[index].mz > cluster_mz + parameters.radius_mz ||
                        peaks[index].mz < cluster_mz - parameters.radius_mz ||
                        peaks[index].rt > cluster_rt + parameters.radius_rt ||
                        peaks[index].rt < cluster_rt - parameters.radius_rt) {
                        peaks[index].cluster_id = -1;
                        metapeak_indexes.erase(metapeak_indexes.begin() + k);
                        calculate_cluster_pos();
                    }
                }
            }
        }

        // Check if this cluster contains the
        // necessary number of peaks per class. If the
        // fraction is not big enough, free the peaks
        // for future clustering.
        // std::vector<size_t> class_ids;
        // std::vector<size_t> class_map;
        // for (const auto& index : metapeak_indexes) {
        // auto& peak = peaks[index];
        // if (peak.cluster_id == cluster_id) {
        // bool id_present = false;
        // for (size_t k = 0; k < class_ids.size(); ++k) {
        // const auto& id = class_ids[k];
        // if (id == peak.class_id) {
        // id_present = true;
        // class_map[k] += 1;
        // break;
        //}
        //}
        // if (!id_present) {
        // class_ids.push_back(peak.class_id);
        // class_map.push_back(1);
        //}
        //}
        //}
        // bool fraction_achieved = false;
        // for (const auto& n_hits : class_map) {
        //// FIXME: Should be n_files_per_class!
        // if ((double)n_hits / (double)parameters.n_files >
        // parameters.fraction) {
        // fraction_achieved = true;
        // break;
        //}
        //}
        std::vector<size_t> class_map(classes.size());
        // ...
        for (const auto& index : metapeak_indexes) {
            const auto& peak = peaks[index];
            for (size_t k = 0; k < classes.size(); ++k) {
                if (peak.class_id == classes[k].id) {
                    class_map[k] += 1;
                    break;
                }
            }
        }
        bool fraction_achieved = false;
        for (size_t k = 0; k < class_map.size(); ++k) {
            const auto& n_hits = class_map[k];
            const auto& n_files = classes[k].file_ids.size();
            if ((double)n_hits / (double)n_files > parameters.fraction) {
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
    std::cout << "SUM PEAKS VISITED: " << avg_peaks_per_iter << std::endl;
    std::cout << "AVG PEAKS VISITED PER ITER: "
              << avg_peaks_per_iter / (double)peaks.size() << std::endl;

    // DEBUG
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
