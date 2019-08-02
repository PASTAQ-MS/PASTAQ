#include "metamatch/metamatch_files.hpp"

bool MetaMatch::Files::Csv::write_clusters(
    std::ostream& stream, const std::vector<MetaMatch::Cluster>& clusters,
    size_t n_files) {
    char cell_delimiter = ' ';
    char line_delimiter = '\n';

    // Prepare the CSV header.
    std::vector<std::string> header_columns = {
        "metapeak",
        "mz",
        "rt",
    };
    for (size_t i = 0; i < n_files; ++i) {
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
        // TODO(alex): Other stats here...
        // file heights
        for (size_t i = 0; i < n_files; ++i) {
            stream << cluster.file_heights[i];
            if (i == n_files - 1) {
                stream << line_delimiter;
            } else {
                stream << cell_delimiter;
            }
        }
    }
    return stream.good();
}

bool MetaMatch::Files::Csv::write_peaks(
    std::ostream& stream, const std::vector<MetaMatch::Peak>& peaks,
    bool include_mpid) {
    char cell_delimiter = ' ';
    char line_delimiter = '\n';

    // Write the CSV header.
    // TODO(alex): Must use only the values that we need. Using the old format
    // for now to test if centroid is working properly.
    std::vector<std::string> header_columns = {
        "N",         "X",        "Y",          "Height",  "Volume",
        "VCentroid", "XSigma",   "YSigma",     "Count",   "LocalBkgnd",
        "SNVolume",  "SNHeight", "SNCentroid", "file_id", "class_id"};
    if (include_mpid) {
        header_columns.push_back("mpid");
    }
    for (size_t i = 0; i < header_columns.size(); ++i) {
        stream << header_columns[i];
        if (i == header_columns.size() - 1) {
            stream << line_delimiter;
        } else {
            stream << cell_delimiter;
        }
    }
    stream.precision(8);

    // Write each peak as a csv row.
    for (size_t i = 0; i < peaks.size(); ++i) {
        auto peak = peaks[i];
        // N
        stream << i
               << cell_delimiter
               // X
               << peak.local_max_mz
               << cell_delimiter
               // Y
               << peak.local_max_rt
               << cell_delimiter
               // Height
               << peak.local_max_height
               << cell_delimiter
               // Volume
               << peak.raw_roi_total_intensity
               << cell_delimiter
               // VCentroid
               << peak.raw_roi_total_intensity
               << cell_delimiter
               // XSigma
               << peak.raw_roi_sigma_mz
               << cell_delimiter
               // YSigma
               << peak.raw_roi_sigma_rt
               << cell_delimiter
               // Count
               << 0
               << cell_delimiter
               // LocalBkgnd
               << peak.slope_descent_border_background
               << cell_delimiter
               // SNVolume
               << (peak.raw_roi_total_intensity /
                   peak.slope_descent_border_background)
               << cell_delimiter
               // SNHeight
               << (peak.local_max_height / peak.slope_descent_border_background)
               << cell_delimiter
               // SNHeight
               << (peak.raw_roi_total_intensity /
                   peak.slope_descent_border_background)
               << cell_delimiter
               // file_id
               << peak.file_id
               << cell_delimiter
               // class_id
               << peak.class_id;
        if (include_mpid) {
            // mpid
            stream << cell_delimiter << peak.cluster_id;
        }
        if (i != peaks.size() - 1) {
            stream << line_delimiter;
        }
    }
    return stream.good();
}
