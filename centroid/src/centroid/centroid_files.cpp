#include "centroid/centroid_files.hpp"
#include "centroid/centroid_serialize.hpp"
#include "grid/grid_serialize.hpp"
#include "utils/serialization.hpp"

bool Centroid::Files::Bpks::write_header(std::ostream &stream,
                                         const Header &header) {
    Serialization::write_uint8(stream, header.spec_version);
    Serialization::write_uint64(stream, header.num_peaks);
    Grid::Serialize::write_parameters(stream, header.grid_params);
    return stream.good();
}

bool Centroid::Files::Bpks::read_header(std::istream &stream, Header *header) {
    Serialization::read_uint8(stream, &header->spec_version);
    Serialization::read_uint64(stream, &header->num_peaks);
    Grid::Serialize::read_parameters(stream, &header->grid_params);
    return stream.good();
}

bool Centroid::Files::Bpks::read_peaks(std::istream &stream,
                                       Grid::Parameters *grid_parameters,
                                       std::vector<Centroid::Peak> *peaks) {
    Header header = {};
    if (!read_header(stream, &header)) {
        return false;
    }
    *grid_parameters = header.grid_params;
    peaks->resize(header.num_peaks);
    for (auto &peak : *peaks) {
        if (!Centroid::Serialize::read_peak(stream, &peak)) {
            return false;
        }
    }
    return stream.good();
}

bool Centroid::Files::Bpks::write_peaks(
    std::ostream &stream, const Grid::Parameters &grid_parameters,
    const std::vector<Centroid::Peak> &peaks) {
    if (!write_header(stream, {0, peaks.size(), grid_parameters})) {
        return false;
    }
    for (const auto &peak : peaks) {
        if (!Centroid::Serialize::write_peak(stream, peak)) {
            return false;
        }
    }
    return stream.good();
}

bool Centroid::Files::Csv::write_peaks(
    std::ostream &stream, const std::vector<Centroid::Peak> &peaks) {
    char cell_delimiter = ' ';
    char line_delimiter = '\n';

    // Write the CSV header.
    // TODO(alex): Must use only the values that we need. Using the old format
    // for now to test if centroid is working properly.
    std::vector<std::string> header_columns = {
        "N",         "X",        "Y",         "Height", "Volume",
        "VCentroid", "XSigma",   "YSigma",    "Count",  "LocalBkgnd",
        "SNVolume",  "SNHeight", "SNCentroid"};
    for (size_t i = 0; i < header_columns.size(); ++i) {
        stream << header_columns[i];
        if (i == header_columns.size() - 1) {
            stream << line_delimiter;
        } else {
            stream << cell_delimiter;
        }
    }

    // Write each peak as a csv row.
    for (size_t i = 0; i < peaks.size(); ++i) {
        auto peak = peaks[i];
        stream << i << cell_delimiter                              // N
               << peak.mz << cell_delimiter                        // X
               << peak.rt << cell_delimiter                        // Y
               << peak.height << cell_delimiter                    // Height
               << peak.total_intensity << cell_delimiter           // Volume
               << peak.total_intensity_centroid << cell_delimiter  // VCentroid
               << peak.sigma_mz << cell_delimiter                  // XSigma
               << peak.sigma_rt << cell_delimiter                  // YSigma
               << peak.points.size() << cell_delimiter             // Count
               << peak.border_background << cell_delimiter         // LocalBkgnd
               << (peak.total_intensity / peak.border_background)
               << cell_delimiter  // SNVolume
               << (peak.height / peak.border_background)
               << cell_delimiter  // SNHeight
               << (peak.total_intensity_centroid /
                   peak.border_background);  // SNHeight
        if (i != peaks.size() - 1) {
            stream << line_delimiter;
        }
    }
    return stream.good();
}
