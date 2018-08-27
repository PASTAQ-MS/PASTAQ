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
    stream.precision(8);

    // Write each peak as a csv row.
    for (size_t i = 0; i < peaks.size(); ++i) {
        auto peak = peaks[i];
        // N
        stream << i
               << cell_delimiter
               // X
               << peak.mz
               << cell_delimiter
               // Y
               << peak.rt
               << cell_delimiter
               // Height
               << peak.height_centroid
               << cell_delimiter
               // Volume
               << peak.total_intensity
               << cell_delimiter
               // VCentroid
               << peak.total_intensity_centroid
               << cell_delimiter
               // XSigma
               << peak.sigma_mz
               << cell_delimiter
               // YSigma
               << peak.sigma_rt
               << cell_delimiter
               // Count
               << peak.points.size()
               << cell_delimiter
               // LocalBkgnd
               << peak.border_background
               << cell_delimiter
               // SNVolume
               << (peak.total_intensity / peak.border_background)
               << cell_delimiter
               // SNHeight
               << (peak.height / peak.border_background)
               << cell_delimiter
               // SNHeight
               << (peak.total_intensity_centroid / peak.border_background);
        if (i != peaks.size() - 1) {
            stream << line_delimiter;
        }
    }
    return stream.good();
}

bool Centroid::Files::Csv::read_peaks(std::istream &stream,
                                      std::vector<Centroid::Peak> *peaks) {
    char cell_delimiter = ' ';
    char line_delimiter = '\n';

    // Read the CSV header.
    // TODO(alex): Must use only the values that we need. Using the old format
    // for now to test if centroid is working properly.
    std::vector<std::string> header_columns = {
        "N",         "X",        "Y",         "Height", "Volume",
        "VCentroid", "XSigma",   "YSigma",    "Count",  "LocalBkgnd",
        "SNVolume",  "SNHeight", "SNCentroid"};
    std::string line;
    if (!std::getline(stream, line, line_delimiter)) {
        return false;
    }
    // Verify that the read header matches the header_columns.
    std::string token;
    std::stringstream token_stream(line);
    int i = 0;
    while (std::getline(token_stream, token, cell_delimiter)) {
        if (token != header_columns[i]) {
            return false;
        };
        ++i;
    }

    // Read the rest of the file and build the list of peaks.
    // TODO(alex): Should we make the parsing independent of the order of the
    // columns?
    while (std::getline(stream, line, line_delimiter)) {
        Centroid::Peak peak = {};
        token_stream = std::stringstream(line);
        // N (Skip).
        std::getline(token_stream, token, cell_delimiter);
        // X
        std::getline(token_stream, token, cell_delimiter);
        if (!(std::istringstream(token) >> peak.mz).eof()) {
            return false;
        };
        peak.mz_centroid = peak.mz;
        // Y
        std::getline(token_stream, token, cell_delimiter);
        if (!(std::istringstream(token) >> peak.rt).eof()) {
            return false;
        };
        peak.rt_centroid = peak.rt;
        // Height
        std::getline(token_stream, token, cell_delimiter);
        if (!(std::istringstream(token) >> peak.height_centroid).eof()) {
            return false;
        };
        peak.height = peak.height_centroid;
        // Volume
        std::getline(token_stream, token, cell_delimiter);
        if (!(std::istringstream(token) >> peak.total_intensity).eof()) {
            return false;
        };
        // VCentroid
        std::getline(token_stream, token, cell_delimiter);
        if (!(std::istringstream(token) >> peak.total_intensity_centroid)
                 .eof()) {
            return false;
        };
        // XSigma
        std::getline(token_stream, token, cell_delimiter);
        if (!(std::istringstream(token) >> peak.sigma_mz).eof()) {
            return false;
        };
        // YSigma
        std::getline(token_stream, token, cell_delimiter);
        if (!(std::istringstream(token) >> peak.sigma_rt).eof()) {
            return false;
        };
        // Count (Skip).
        std::getline(token_stream, token, cell_delimiter);
        // LocalBkgnd
        std::getline(token_stream, token, cell_delimiter);
        if (!(std::istringstream(token) >> peak.border_background).eof()) {
            return false;
        };
        // SNVolume (Skip)
        // SNHeight (Skip)
        // SNHeight (Skip)

        // Add the peak to the list.
        peaks->push_back(peak);
    }

    return stream.good() || stream.eof();
}
