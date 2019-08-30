#include <sstream>

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

bool Centroid::Files::Bpks::read_peaks(std::istream &stream,
                                       std::vector<Centroid::Peak> *peaks) {
    uint64_t num_peaks = 0;
    Serialization::read_uint64(stream, &num_peaks);
    peaks->resize(num_peaks);
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

bool Centroid::Files::Bpks::write_peaks(
    std::ostream &stream, const std::vector<Centroid::Peak> &peaks) {
    Serialization::write_uint64(stream, peaks.size());
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
               // TODO: RETURN raw_roi_total_intensity here? other?
               //<< peak.raw_roi_total_intensity
               << peak.local_max_height
               << cell_delimiter
               // Volume
               << peak.raw_roi_total_intensity
               << cell_delimiter
               // VCentroid
               << 0
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
               //<< peak.slope_descent_border_background
               << 0
               << cell_delimiter
               // SNVolume
               << 0
               << cell_delimiter
               // SNHeight
               << 0
               << cell_delimiter
               // SNCentroid
               << 0;
        if (i != peaks.size() - 1) {
            stream << line_delimiter;
        }
    }
    return stream.good();
}

bool Centroid::Files::Csv::read_peaks(std::istream &stream,
                                      std::vector<Centroid::Peak> *peaks) {
    return stream.good() || stream.eof();
}
