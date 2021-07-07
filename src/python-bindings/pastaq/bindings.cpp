#include <cassert>
#include <fstream>
#include <iostream>
#include <thread>
#include <tuple>

#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "centroid/centroid.hpp"
#include "centroid/centroid_serialize.hpp"
#include "feature_detection/feature_detection.hpp"
#include "feature_detection/feature_detection_serialize.hpp"
#include "grid/grid.hpp"
#include "grid/grid_serialize.hpp"
#include "link/link.hpp"
#include "link/link_serialize.hpp"
#include "metamatch/metamatch.hpp"
#include "metamatch/metamatch_serialize.hpp"
#include "protein_inference/protein_inference.hpp"
#include "protein_inference/protein_inference_serialize.hpp"
#include "raw_data/raw_data.hpp"
#include "raw_data/raw_data_serialize.hpp"
#include "raw_data/xml_reader.hpp"
#include "utils/compression.hpp"
#include "utils/search.hpp"
#include "utils/serialization.hpp"
#include "warp2d/warp2d.hpp"
#include "warp2d/warp2d_serialize.hpp"

namespace py = pybind11;

namespace PythonAPI {
RawData::RawData read_mzxml(std::string &input_file, double min_mz,
                            double max_mz, double min_rt, double max_rt,
                            std::string instrument_type_str,
                            double resolution_ms1, double resolution_msn,
                            double reference_mz, double fwhm_rt,
                            std::string polarity_str, size_t ms_level) {
    pybind11::gil_scoped_release release;
    // Setup infinite range if no point was specified.
    min_rt = min_rt < 0 ? 0 : min_rt;
    max_rt = max_rt < 0 ? std::numeric_limits<double>::infinity() : max_rt;
    min_mz = min_mz < 0 ? 0 : min_mz;
    max_mz = max_mz < 0 ? std::numeric_limits<double>::infinity() : max_mz;

    // Parse the instrument type.
    auto instrument_type = Instrument::UNKNOWN;
    for (auto &ch : instrument_type_str) {
        ch = std::tolower(ch);
    }
    if (instrument_type_str == "orbitrap") {
        instrument_type = Instrument::ORBITRAP;
    } else if (instrument_type_str == "tof") {
        instrument_type = Instrument::TOF;
    } else if (instrument_type_str == "quad" || instrument_type_str == "quadrupole") {
        instrument_type = Instrument::QUAD;
    } else if (instrument_type_str == "fticr" || instrument_type_str == "ft-icr" ) {
        instrument_type = Instrument::FTICR;
    } else {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "the given instrument is not supported";
        throw std::invalid_argument(error_stream.str());
    }
    // Parse the polarity.
    auto polarity = Polarity::BOTH;
    for (auto &ch : polarity_str) {
        ch = std::tolower(ch);
    }
    if (polarity_str == "" || polarity_str == "both" || polarity_str == "+-" ||
        polarity_str == "-+") {
        polarity = Polarity::BOTH;
    } else if (polarity_str == "+" || polarity_str == "pos" ||
               polarity_str == "positive") {
        polarity = Polarity::POSITIVE;
    } else if (polarity_str == "-" || polarity_str == "neg" ||
               polarity_str == "negative") {
        polarity = Polarity::NEGATIVE;
    } else {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "the given polarity is not supported. choose "
                        "between '+', '-', 'both' (default)";
        throw std::invalid_argument(error_stream.str());
    }

    // Sanity check the min/max rt/mz.
    if (min_rt >= max_rt) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: min_rt >= max_rt (min_rt: " << min_rt
                     << ", max_rt: " << max_rt << ")";
        throw std::invalid_argument(error_stream.str());
    }
    if (min_mz >= max_mz) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: min_mz >= max_mz (min_mz: " << min_mz
                     << ", max_mz: " << max_mz << ")";
        throw std::invalid_argument(error_stream.str());
    }

    // Open file stream.
    std::ifstream stream;
    stream.open(input_file);
    if (!stream) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }

    auto raw_data = XmlReader::read_mzxml(
        stream, min_mz, max_mz, min_rt, max_rt, instrument_type, resolution_ms1,
        resolution_msn, reference_mz, polarity, ms_level);
    if (!raw_data) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: an error occurred when reading the file"
                     << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    raw_data.value().fwhm_rt = fwhm_rt;
    pybind11::gil_scoped_acquire acquire;

    return raw_data.value();
}

RawData::RawData read_mzml(std::string &input_file, double min_mz,
                           double max_mz, double min_rt, double max_rt,
                           std::string instrument_type_str,
                           double resolution_ms1, double resolution_msn,
                           double reference_mz, double fwhm_rt,
                           std::string polarity_str, size_t ms_level) {
    pybind11::gil_scoped_release release;
    // Setup infinite range if no point was specified.
    min_rt = min_rt < 0 ? 0 : min_rt;
    max_rt = max_rt < 0 ? std::numeric_limits<double>::infinity() : max_rt;
    min_mz = min_mz < 0 ? 0 : min_mz;
    max_mz = max_mz < 0 ? std::numeric_limits<double>::infinity() : max_mz;

    // Parse the instrument type.
    auto instrument_type = Instrument::UNKNOWN;
    for (auto &ch : instrument_type_str) {
        ch = std::tolower(ch);
    }
    if (instrument_type_str == "orbitrap") {
        instrument_type = Instrument::ORBITRAP;
    } else if (instrument_type_str == "tof") {
        instrument_type = Instrument::TOF;
    } else if (instrument_type_str == "quad") {
        instrument_type = Instrument::QUAD;
    } else if (instrument_type_str == "fticr") {
        instrument_type = Instrument::FTICR;
    } else {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "the given instrument is not supported";
        throw std::invalid_argument(error_stream.str());
    }
    // Parse the polarity.
    auto polarity = Polarity::BOTH;
    for (auto &ch : polarity_str) {
        ch = std::tolower(ch);
    }
    if (polarity_str == "" || polarity_str == "both" || polarity_str == "+-" ||
        polarity_str == "-+") {
        polarity = Polarity::BOTH;
    } else if (polarity_str == "+" || polarity_str == "pos" ||
               polarity_str == "positive") {
        polarity = Polarity::POSITIVE;
    } else if (polarity_str == "-" || polarity_str == "neg" ||
               polarity_str == "negative") {
        polarity = Polarity::NEGATIVE;
    } else {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "the given polarity is not supported. choose "
                        "between '+', '-', 'both' (default)";
        throw std::invalid_argument(error_stream.str());
    }

    // Sanity check the min/max rt/mz.
    if (min_rt >= max_rt) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: min_rt >= max_rt (min_rt: " << min_rt
                     << ", max_rt: " << max_rt << ")";
        throw std::invalid_argument(error_stream.str());
    }
    if (min_mz >= max_mz) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: min_mz >= max_mz (min_mz: " << min_mz
                     << ", max_mz: " << max_mz << ")";
        throw std::invalid_argument(error_stream.str());
    }

    // Open file stream.
    std::ifstream stream;
    stream.open(input_file);
    if (!stream) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }

    auto raw_data = XmlReader::read_mzml(
        stream, min_mz, max_mz, min_rt, max_rt, instrument_type, resolution_ms1,
        resolution_msn, reference_mz, polarity, ms_level);
    if (!raw_data) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: an error occurred when reading the file"
                     << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    raw_data.value().fwhm_rt = fwhm_rt;

    pybind11::gil_scoped_acquire acquire;
    return raw_data.value();
}

Xic::Xic xic(const RawData::RawData &raw_data, double min_mz, double max_mz,
             double min_rt, double max_rt, std::string method_str) {
    pybind11::gil_scoped_release release;
    // Parse the instrument type.
    auto method = Xic::UNKNOWN;
    for (auto &ch : method_str) {
        ch = std::tolower(ch);
    }
    if (method_str == "max") {
        method = Xic::MAX;
    } else if (method_str == "sum") {
        method = Xic::SUM;
    } else {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "the given xic method is not supported";
        throw std::invalid_argument(error_stream.str());
    }
    pybind11::gil_scoped_acquire acquire;
    return RawData::xic(raw_data, min_mz, max_mz, min_rt, max_rt, method);
}

Grid::Grid resample(const RawData::RawData &raw_data, uint64_t num_samples_mz,
                    uint64_t num_samples_rt, double smoothing_coef_mz,
                    double smoothing_coef_rt) {
    pybind11::gil_scoped_release release;
    auto params = Grid::ResampleParams{};
    params.num_samples_mz = num_samples_mz;
    params.num_samples_rt = num_samples_rt;
    params.smoothing_coef_mz = smoothing_coef_mz;
    params.smoothing_coef_rt = smoothing_coef_rt;
    auto grid =  Grid::resample(raw_data, params);
    pybind11::gil_scoped_acquire acquire;
    return grid;
}

std::string to_string(const Instrument::Type &instrument_type) {
    switch (instrument_type) {
        case Instrument::QUAD:
            return "QUAD";
        case Instrument::TOF:
            return "TOF";
        case Instrument::FTICR:
            return "FTICR";
        case Instrument::ORBITRAP:
            return "ORBITRAP";
        default:
            return "UNKNOWN";
    };
}

std::string to_string(const Polarity::Type &polarity) {
    switch (polarity) {
        case Polarity::POSITIVE:
            return "POSITIVE";
        case Polarity::NEGATIVE:
            return "NEGATIVE";
        case Polarity::BOTH:
            return "BOTH";
        default:
            return "UNKNOWN";
    };
}

std::string to_string(const Xic::Method &method) {
    switch (method) {
        case Xic::Method::MAX:
            return "MAX";
        case Xic::Method::SUM:
            return "SUM";
        default:
            return "UNKNOWN";
    };
}

Warp2D::TimeMap calculate_time_map(
    const std::vector<Centroid::Peak> &ref_peaks,
    const std::vector<Centroid::Peak> &source_peaks, int64_t slack,
    int64_t window_size, int64_t num_points, double rt_expand_factor,
    int64_t peaks_per_window) {
    pybind11::gil_scoped_release release;
    // TODO(alex): Validate the parameters and throw an error if
    // appropriate.
    Warp2D::Parameters parameters = {slack, window_size, num_points,
                                     peaks_per_window, rt_expand_factor};
    auto time_map =
        Warp2D::calculate_time_map(ref_peaks, source_peaks, parameters,
                                   std::thread::hardware_concurrency());
    pybind11::gil_scoped_acquire acquire;
    return time_map;
}

// TODO: Where should this function go?
struct SimilarityResults {
    double self_a;
    double self_b;
    double overlap;
    double geometric_ratio;
    double mean_ratio;
};
SimilarityResults find_similarity(std::vector<Centroid::Peak> &peak_list_a,
                                  std::vector<Centroid::Peak> &peak_list_b,
                                  size_t n_peaks) {
    pybind11::gil_scoped_release release;
    auto sort_peaks = [](const Centroid::Peak &p1,
                         const Centroid::Peak &p2) -> bool {
        return (p2.fitted_height < p1.fitted_height);
    };
    std::sort(peak_list_a.begin(), peak_list_a.end(), sort_peaks);
    std::sort(peak_list_b.begin(), peak_list_b.end(), sort_peaks);
    if (peak_list_a.size() > n_peaks) {
        peak_list_a.resize(n_peaks);
    }
    if (peak_list_b.size() > n_peaks) {
        peak_list_b.resize(n_peaks);
    }
    SimilarityResults results = {};
    results.self_a = Centroid::cumulative_overlap(peak_list_a, peak_list_a);
    results.self_b = Centroid::cumulative_overlap(peak_list_b, peak_list_b);
    results.overlap = Centroid::cumulative_overlap(peak_list_a, peak_list_b);
    results.geometric_ratio = 0;
    results.mean_ratio = 0;
    if (results.self_a != 0 && results.self_b != 0) {
        // Overlap / (GeometricMean(self_a, self_b))
        results.geometric_ratio =
            results.overlap / std::sqrt(results.self_a * results.self_b);
        // Harmonic mean of the ratios between
        // self_similarity/overlap_similarity
        results.mean_ratio =
            2 * results.overlap / (results.self_a + results.self_b);
    }
    pybind11::gil_scoped_acquire acquire;
    return results;
}

void write_raw_data(const RawData::RawData &raw_data,
                    std::string &output_file) {
    pybind11::gil_scoped_release release;
    // Open file stream.
    Compression::DeflateStream stream;
    stream.open(output_file);
    if (!stream) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't open output file" << output_file;
        throw std::invalid_argument(error_stream.str());
    }

    if (!RawData::Serialize::write_raw_data(stream, raw_data)) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream
            << "error: couldn't write the raw_data into the output file"
            << output_file;
        throw std::invalid_argument(error_stream.str());
    }
    pybind11::gil_scoped_acquire acquire;
}

RawData::RawData read_raw_data(std::string &input_file) {
    pybind11::gil_scoped_release release;
    // Open file stream.
    Compression::InflateStream stream;
    stream.open(input_file);
    if (!stream) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }

    RawData::RawData raw_data;
    if (!RawData::Serialize::read_raw_data(stream, &raw_data)) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't write the raw_data into the input file"
                     << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    pybind11::gil_scoped_acquire acquire;
    return raw_data;
}

void write_grid(const Grid::Grid &grid, std::string &output_file) {
    pybind11::gil_scoped_release release;
    // Open file stream.
    Compression::DeflateStream stream;
    stream.open(output_file);
    if (!stream) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't open output file" << output_file;
        throw std::invalid_argument(error_stream.str());
    }

    if (!Grid::Serialize::write_grid(stream, grid)) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't write the grid into the output file"
                     << output_file;
        throw std::invalid_argument(error_stream.str());
    }
    pybind11::gil_scoped_acquire acquire;
}

Grid::Grid read_grid(std::string &input_file) {
    pybind11::gil_scoped_release release;
    // Open file stream.
    Compression::InflateStream stream;
    stream.open(input_file);
    if (!stream) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }

    Grid::Grid grid;
    if (!Grid::Serialize::read_grid(stream, &grid)) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't write the grid into the input file"
                     << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    pybind11::gil_scoped_acquire acquire;
    return grid;
}

void write_time_map(const Warp2D::TimeMap &time_map, std::string &output_file) {
    pybind11::gil_scoped_release release;
    // Open file stream.
    Compression::DeflateStream stream;
    stream.open(output_file);
    if (!stream) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't open output file" << output_file;
        throw std::invalid_argument(error_stream.str());
    }

    if (!Warp2D::Serialize::write_time_map(stream, time_map)) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream
            << "error: couldn't write the time_map into the output file"
            << output_file;
        throw std::invalid_argument(error_stream.str());
    }
    pybind11::gil_scoped_acquire acquire;
}

Warp2D::TimeMap read_time_map(std::string &input_file) {
    pybind11::gil_scoped_release release;
    // Open file stream.
    Compression::InflateStream stream;
    stream.open(input_file);
    if (!stream) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }

    Warp2D::TimeMap time_map;
    if (!Warp2D::Serialize::read_time_map(stream, &time_map)) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't write the time_map into the input file"
                     << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    pybind11::gil_scoped_acquire acquire;
    return time_map;
}

void write_peaks(const std::vector<Centroid::Peak> &peaks,
                 std::string &output_file) {
    pybind11::gil_scoped_release release;
    // Open file stream.
    Compression::DeflateStream stream;
    stream.open(output_file);
    if (!stream) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't open output file" << output_file;
        throw std::invalid_argument(error_stream.str());
    }

    if (!Centroid::Serialize::write_peaks(stream, peaks)) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't write the peaks into the output file"
                     << output_file;
        throw std::invalid_argument(error_stream.str());
    }
    pybind11::gil_scoped_acquire acquire;
}

std::vector<Centroid::Peak> read_peaks(std::string &input_file) {
    pybind11::gil_scoped_release release;
    // Open file stream.
    Compression::InflateStream stream;
    stream.open(input_file);
    if (!stream) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }

    std::vector<Centroid::Peak> peaks;
    if (!Centroid::Serialize::read_peaks(stream, &peaks)) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't write the peaks into the input file"
                     << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    pybind11::gil_scoped_acquire acquire;
    return peaks;
}

std::vector<ProteinInference::InferredProtein> read_inferred_proteins(
    std::string &input_file) {
    pybind11::gil_scoped_release release;
    // Open file stream.
    Compression::InflateStream stream;
    stream.open(input_file);
    if (!stream) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }

    std::vector<ProteinInference::InferredProtein> inferred_proteins;
    if (!ProteinInference::Serialize::read_inferred_proteins(
            stream, &inferred_proteins)) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream
            << "error: couldn't write the inferred_proteins into the input file"
            << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    pybind11::gil_scoped_acquire acquire;
    return inferred_proteins;
}

void write_inferred_proteins(
    const std::vector<ProteinInference::InferredProtein> &inferred_proteins,
    std::string &output_file) {
    pybind11::gil_scoped_release release;
    // Open file stream.
    Compression::DeflateStream stream;
    stream.open(output_file);
    if (!stream) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't open output file" << output_file;
        throw std::invalid_argument(error_stream.str());
    }

    if (!ProteinInference::Serialize::write_inferred_proteins(
            stream, inferred_proteins)) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't write the inferred_proteins into the "
                        "output file"
                     << output_file;
        throw std::invalid_argument(error_stream.str());
    }
    pybind11::gil_scoped_acquire acquire;
}

void write_feature_clusters(
    const std::vector<MetaMatch::FeatureCluster> &feature_clusters,
    std::string &output_file) {
    pybind11::gil_scoped_release release;
    // Open file stream.
    Compression::DeflateStream stream;
    stream.open(output_file);
    if (!stream) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't open output file" << output_file;
        throw std::invalid_argument(error_stream.str());
    }

    if (!MetaMatch::Serialize::write_feature_clusters(stream,
                                                      feature_clusters)) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream
            << "error: couldn't write the feature_clusters into the output file"
            << output_file;
        throw std::invalid_argument(error_stream.str());
    }
    pybind11::gil_scoped_acquire acquire;
}

std::vector<MetaMatch::FeatureCluster> read_feature_clusters(
    std::string &input_file) {
    pybind11::gil_scoped_release release;
    // Open file stream.
    Compression::InflateStream stream;
    stream.open(input_file);
    if (!stream) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }

    std::vector<MetaMatch::FeatureCluster> feature_clusters;
    if (!MetaMatch::Serialize::read_feature_clusters(stream,
                                                     &feature_clusters)) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream
            << "error: couldn't write the feature_clusters into the input file"
            << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    pybind11::gil_scoped_acquire acquire;
    return feature_clusters;
}

void write_peak_clusters(
    const std::vector<MetaMatch::PeakCluster> &peak_clusters,
    std::string &output_file) {
    pybind11::gil_scoped_release release;
    // Open file stream.
    Compression::DeflateStream stream;
    stream.open(output_file);
    if (!stream) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't open output file" << output_file;
        throw std::invalid_argument(error_stream.str());
    }

    if (!MetaMatch::Serialize::write_peak_clusters(stream, peak_clusters)) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream
            << "error: couldn't write the peak_clusters into the output file"
            << output_file;
        throw std::invalid_argument(error_stream.str());
    }
    pybind11::gil_scoped_acquire acquire;
}

std::vector<MetaMatch::PeakCluster> read_peak_clusters(
    std::string &input_file) {
    pybind11::gil_scoped_release release;
    // Open file stream.
    Compression::InflateStream stream;
    stream.open(input_file);
    if (!stream) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }

    std::vector<MetaMatch::PeakCluster> peak_clusters;
    if (!MetaMatch::Serialize::read_peak_clusters(stream, &peak_clusters)) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream
            << "error: couldn't write the peak_clusters into the input file"
            << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    pybind11::gil_scoped_acquire acquire;
    return peak_clusters;
}

void write_features(const std::vector<FeatureDetection::Feature> &features,
                    std::string &output_file) {
    pybind11::gil_scoped_release release;
    // Open file stream.
    Compression::DeflateStream stream;
    stream.open(output_file);
    if (!stream) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't open output file" << output_file;
        throw std::invalid_argument(error_stream.str());
    }

    if (!FeatureDetection::Serialize::write_features(stream, features)) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream
            << "error: couldn't write the features into the output file"
            << output_file;
        throw std::invalid_argument(error_stream.str());
    }
    pybind11::gil_scoped_acquire acquire;
}

std::vector<FeatureDetection::Feature> read_features(std::string &input_file) {
    pybind11::gil_scoped_release release;
    // Open file stream.
    Compression::InflateStream stream;
    stream.open(input_file);
    if (!stream) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }

    std::vector<FeatureDetection::Feature> features;
    if (!FeatureDetection::Serialize::read_features(stream, &features)) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't write the features into the input file"
                     << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    pybind11::gil_scoped_acquire acquire;
    return features;
}

void write_ident_data(const IdentData::IdentData &ident_data,
                      std::string &output_file) {
    pybind11::gil_scoped_release release;
    // Open file stream.
    Compression::DeflateStream stream;
    stream.open(output_file);
    if (!stream) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't open output file" << output_file;
        throw std::invalid_argument(error_stream.str());
    }

    if (!IdentData::Serialize::write_ident_data(stream, ident_data)) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream
            << "error: couldn't write the ident_data into the output file"
            << output_file;
        throw std::invalid_argument(error_stream.str());
    }
    pybind11::gil_scoped_acquire acquire;
}

IdentData::IdentData read_ident_data(std::string &input_file) {
    pybind11::gil_scoped_release release;
    // Open file stream.
    Compression::InflateStream stream;
    stream.open(input_file);
    if (!stream) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }

    IdentData::IdentData ident_data;
    if (!IdentData::Serialize::read_ident_data(stream, &ident_data)) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream
            << "error: couldn't write the ident_data into the input file"
            << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    pybind11::gil_scoped_acquire acquire;
    return ident_data;
}

void write_linked_msms(const std::vector<Link::LinkedMsms> &linked_msms,
                       std::string &output_file) {
    pybind11::gil_scoped_release release;
    // Open file stream.
    Compression::DeflateStream stream;
    stream.open(output_file);
    if (!stream) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't open output file" << output_file;
        throw std::invalid_argument(error_stream.str());
    }

    if (!Link::Serialize::write_linked_msms_table(stream, linked_msms)) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream
            << "error: couldn't write the linked_msms into the output file"
            << output_file;
        throw std::invalid_argument(error_stream.str());
    }
    pybind11::gil_scoped_acquire acquire;
}

std::vector<Link::LinkedMsms> read_linked_msms(std::string &input_file) {
    pybind11::gil_scoped_release release;
    // Open file stream.
    Compression::InflateStream stream;
    stream.open(input_file);
    if (!stream) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }

    std::vector<Link::LinkedMsms> linked_msms;
    if (!Link::Serialize::read_linked_msms_table(stream, &linked_msms)) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream
            << "error: couldn't write the linked_msms into the input file"
            << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    pybind11::gil_scoped_acquire acquire;
    return linked_msms;
}

void write_linked_psm(const std::vector<Link::LinkedPsm> &linked_psm,
                      std::string &output_file) {
    pybind11::gil_scoped_release release;
    // Open file stream.
    Compression::DeflateStream stream;
    stream.open(output_file);
    if (!stream) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't open output file" << output_file;
        throw std::invalid_argument(error_stream.str());
    }

    if (!Link::Serialize::write_linked_psm_table(stream, linked_psm)) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream
            << "error: couldn't write the linked_psm into the output file"
            << output_file;
        throw std::invalid_argument(error_stream.str());
    }
    pybind11::gil_scoped_acquire acquire;
}

std::vector<Link::LinkedPsm> read_linked_psm(std::string &input_file) {
    pybind11::gil_scoped_release release;
    // Open file stream.
    Compression::InflateStream stream;
    stream.open(input_file);
    if (!stream) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }

    std::vector<Link::LinkedPsm> linked_psm;
    if (!Link::Serialize::read_linked_psm_table(stream, &linked_psm)) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream
            << "error: couldn't write the linked_psm into the input file"
            << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    pybind11::gil_scoped_acquire acquire;
    return linked_psm;
}

IdentData::IdentData read_mzidentml(std::string &input_file, bool ignore_decoy,
                                    bool require_threshold,
                                    bool max_rank_only,
                                    double min_mz,
                                    double max_mz,
                                    double min_rt,
                                    double max_rt) {
    pybind11::gil_scoped_release release;
    // Setup infinite range if no point was specified.
    min_rt = min_rt < 0 ? 0 : min_rt;
    max_rt = max_rt < 0 ? std::numeric_limits<double>::infinity() : max_rt;
    min_mz = min_mz < 0 ? 0 : min_mz;
    max_mz = max_mz < 0 ? std::numeric_limits<double>::infinity() : max_mz;

    // Open file stream.
    std::ifstream stream;
    stream.open(input_file);
    if (!stream) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    auto ident_data = XmlReader::read_mzidentml(stream, ignore_decoy, require_threshold,
                                     max_rank_only, min_mz, max_mz, min_rt, max_rt);
    pybind11::gil_scoped_acquire acquire;
    return ident_data;
}

std::vector<MetaMatch::FeatureCluster> find_feature_clusters(
    std::vector<uint64_t> group_ids,
    std::vector<std::vector<FeatureDetection::Feature>> features,
    double keep_perc, double intensity_threshold, double n_sig_mz,
    double n_sig_rt) {
    pybind11::gil_scoped_release release;
    if (group_ids.size() != features.size()) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: the length of groups and features don't match";
        throw std::invalid_argument(error_stream.str());
    }
    auto clusters = MetaMatch::find_feature_clusters(group_ids, features, keep_perc,
                                            intensity_threshold, n_sig_mz,
                                            n_sig_rt);
    pybind11::gil_scoped_acquire acquire;
    return clusters;
}

std::vector<MetaMatch::PeakCluster> find_peak_clusters(
    std::vector<uint64_t> group_ids,
    std::vector<std::vector<Centroid::Peak>> peaks,
    double keep_perc, double intensity_threshold, double n_sig_mz,
    double n_sig_rt) {
    pybind11::gil_scoped_release release;
    if (group_ids.size() != peaks.size()) {
        pybind11::gil_scoped_acquire acquire;
        std::ostringstream error_stream;
        error_stream << "error: the length of groups and peaks don't match";
        throw std::invalid_argument(error_stream.str());
    }
    auto clusters = MetaMatch::find_peak_clusters(group_ids, peaks, keep_perc,
                                            intensity_threshold, n_sig_mz,
                                            n_sig_rt);
    pybind11::gil_scoped_acquire acquire;
    return clusters;
}

}  // namespace PythonAPI

PYBIND11_MODULE(pastaq, m) {
    // Documentation.
    m.doc() = "pastaq documentation";

    // Structs.
    py::class_<RawData::PrecursorInformation>(m, "PrecursorInformation")
        .def_readonly("id", &RawData::PrecursorInformation::scan_number)
        .def_readonly("charge", &RawData::PrecursorInformation::charge)
        .def_readonly("mz", &RawData::PrecursorInformation::mz)
        .def_readonly("intensity", &RawData::PrecursorInformation::intensity)
        .def_readonly("window_wideness",
                      &RawData::PrecursorInformation::window_wideness);

    py::class_<RawData::Scan>(m, "Scan")
        .def_readonly("scan_number", &RawData::Scan::scan_number)
        .def_readonly("ms_level", &RawData::Scan::ms_level)
        .def_readonly("num_points", &RawData::Scan::num_points)
        .def_readonly("retention_time", &RawData::Scan::retention_time)
        .def_readonly("mz", &RawData::Scan::mz)
        .def_readonly("intensity", &RawData::Scan::intensity)
        .def_readonly("polarity", &RawData::Scan::polarity)
        .def_readonly("precursor_information",
                      &RawData::Scan::precursor_information)
        .def("__repr__", [](const RawData::Scan &s) {
            auto msg = "Scan <id: " + std::to_string(s.scan_number) +
                       ", ms_level: " + std::to_string(s.ms_level) +
                       ", retention_time: " + std::to_string(s.retention_time) +
                       ", num_points: " + std::to_string(s.num_points) +
                       ", polarity: " + PythonAPI::to_string(s.polarity);
            if (s.ms_level != 1) {
                msg += ", precursor_id: " +
                       std::to_string(s.precursor_information.scan_number);
                msg += ", precursor_charge: " +
                       std::to_string(s.precursor_information.charge);
                msg += ", precursor_mz: " +
                       std::to_string(s.precursor_information.mz);
                msg += ", precursor_intensity: " +
                       std::to_string(s.precursor_information.intensity);
                msg += ", precursor_window_wideness: " +
                       std::to_string(s.precursor_information.window_wideness);
            }
            return msg + ">";
        });

    py::class_<Instrument::Type>(m, "Instrument")
        .def("__repr__", [](const Instrument::Type &instrument_type) {
            return PythonAPI::to_string(instrument_type);
        });

    py::class_<Polarity::Type>(m, "Polarity")
        .def("__repr__", [](const Polarity::Type &polarity) {
            return PythonAPI::to_string(polarity);
        });

    py::class_<RawData::RawData>(m, "RawData")
        .def_readonly("scans", &RawData::RawData::scans)
        .def_readonly("fwhm_rt", &RawData::RawData::fwhm_rt)
        .def_readonly("instrument_type", &RawData::RawData::instrument_type)
        .def_readonly("resolution_ms1", &RawData::RawData::resolution_ms1)
        .def_readonly("resolution_msn", &RawData::RawData::resolution_msn)
        .def_readonly("reference_mz", &RawData::RawData::reference_mz)
        .def_readonly("min_mz", &RawData::RawData::min_mz)
        .def_readonly("max_mz", &RawData::RawData::max_mz)
        .def_readonly("min_rt", &RawData::RawData::min_rt)
        .def_readonly("max_rt", &RawData::RawData::max_rt)
        .def("theoretical_fwhm", &RawData::theoretical_fwhm, py::arg("mz"))
        .def("dump", &PythonAPI::write_raw_data)
        .def("raw_points", &RawData::raw_points,
             "Get the raw data points on the square region defined by "
             "min/max_mz/rt",
             py::arg("min_mz"), py::arg("max_mz"), py::arg("min_rt"),
             py::arg("max_rt"))
        .def("__repr__", [](const RawData::RawData &rd) {
            return "RawData:\n> instrument_type: " +
                   PythonAPI::to_string(rd.instrument_type) +
                   "\n> resolution_ms1: " + std::to_string(rd.resolution_ms1) +
                   "\n> resolution_msn: " + std::to_string(rd.resolution_msn) +
                   "\n> reference_mz: " + std::to_string(rd.reference_mz) +
                   "\n> min_mz: " + std::to_string(rd.min_mz) +
                   "\n> max_mz: " + std::to_string(rd.max_mz) +
                   "\n> min_rt: " + std::to_string(rd.min_rt) +
                   "\n> max_rt: " + std::to_string(rd.max_rt) +
                   "\n> number of scans: " + std::to_string(rd.scans.size());
        });

    py::class_<Grid::Grid>(m, "Grid")
        .def_readonly("n", &Grid::Grid::n)
        .def_readonly("m", &Grid::Grid::m)
        .def_readonly("data", &Grid::Grid::data)
        .def_readonly("bins_mz", &Grid::Grid::bins_mz)
        .def_readonly("bins_rt", &Grid::Grid::bins_rt)
        .def("dump", &PythonAPI::write_grid)
        .def("subset", &Grid::subset)
        .def("__repr__", [](const Grid::Grid &s) {
            return "Grid <n: " + std::to_string(s.n) +
                   ", m: " + std::to_string(s.m) +
                   ", k: " + std::to_string(s.k) +
                   ", t: " + std::to_string(s.t) +
                   ", min_mz: " + std::to_string(s.min_mz) +
                   ", max_mz: " + std::to_string(s.max_mz) +
                   ", min_rt: " + std::to_string(s.min_rt) +
                   ", max_rt: " + std::to_string(s.max_rt) + ">";
        });

    py::class_<RawData::RawPoints>(m, "RawPoints")
        .def_readonly("rt", &RawData::RawPoints::rt)
        .def_readonly("mz", &RawData::RawPoints::mz)
        .def_readonly("intensity", &RawData::RawPoints::intensity);

    py::class_<Xic::Xic>(m, "Xic")
        .def_readonly("retention_time", &Xic::Xic::retention_time)
        .def_readonly("intensity", &Xic::Xic::intensity)
        .def("__repr__", [](const Xic::Xic &s) {
            return "Xic <method: " + PythonAPI::to_string(s.method) +
                   ", min_mz: " + std::to_string(s.min_mz) +
                   ", max_mz: " + std::to_string(s.max_mz) +
                   ", min_rt: " + std::to_string(s.min_rt) +
                   ", max_rt: " + std::to_string(s.max_rt) + ">";
        });

    py::class_<Centroid::Peak>(m, "Peak")
        .def_readonly("id", &Centroid::Peak::id)
        .def_readonly("local_max_mz", &Centroid::Peak::local_max_mz)
        .def_readonly("local_max_rt", &Centroid::Peak::local_max_rt)
        .def_readonly("local_max_height", &Centroid::Peak::local_max_height)
        .def_readonly("rt_delta", &Centroid::Peak::rt_delta)
        .def_readonly("roi_min_mz", &Centroid::Peak::roi_min_mz)
        .def_readonly("roi_max_mz", &Centroid::Peak::roi_max_mz)
        .def_readonly("roi_min_rt", &Centroid::Peak::roi_min_rt)
        .def_readonly("roi_max_rt", &Centroid::Peak::roi_max_rt)
        .def_readonly("raw_roi_mean_mz", &Centroid::Peak::raw_roi_mean_mz)
        .def_readonly("raw_roi_mean_rt", &Centroid::Peak::raw_roi_mean_rt)
        .def_readonly("raw_roi_sigma_mz", &Centroid::Peak::raw_roi_sigma_mz)
        .def_readonly("raw_roi_sigma_rt", &Centroid::Peak::raw_roi_sigma_rt)
        .def_readonly("raw_roi_skewness_mz",
                      &Centroid::Peak::raw_roi_skewness_mz)
        .def_readonly("raw_roi_skewness_rt",
                      &Centroid::Peak::raw_roi_skewness_rt)
        .def_readonly("raw_roi_kurtosis_mz",
                      &Centroid::Peak::raw_roi_kurtosis_mz)
        .def_readonly("raw_roi_kurtosis_rt",
                      &Centroid::Peak::raw_roi_kurtosis_rt)
        .def_readonly("raw_roi_total_intensity",
                      &Centroid::Peak::raw_roi_total_intensity)
        .def_readonly("raw_roi_max_height", &Centroid::Peak::raw_roi_max_height)
        .def_readonly("raw_roi_num_points", &Centroid::Peak::raw_roi_num_points)
        .def_readonly("raw_roi_num_scans", &Centroid::Peak::raw_roi_num_scans)
        .def_readonly("fitted_height", &Centroid::Peak::fitted_height)
        .def_readonly("fitted_mz", &Centroid::Peak::fitted_mz)
        .def_readonly("fitted_rt", &Centroid::Peak::fitted_rt)
        .def_readonly("fitted_sigma_mz", &Centroid::Peak::fitted_sigma_mz)
        .def_readonly("fitted_sigma_rt", &Centroid::Peak::fitted_sigma_rt)
        .def_readonly("fitted_volume", &Centroid::Peak::fitted_volume)
        .def(
            "xic",
            [](const Centroid::Peak &peak, const RawData::RawData &raw_data,
               std::string method) {
                return PythonAPI::xic(raw_data, peak.roi_min_mz,
                                      peak.roi_max_mz, peak.roi_min_rt,
                                      peak.roi_max_rt, method);
            },
            "Get the raw data points on the square region defined by "
            "min/max_mz/rt",
            py::arg("raw_data"), py::arg("method") = "sum")
        .def("__repr__", [](const Centroid::Peak &p) {
            std::string ret = "";
            ret += "Peak <id: " + std::to_string(p.id);
            ret += ", local_max_mz: " + std::to_string(p.local_max_mz);
            ret += ", local_max_rt: " + std::to_string(p.local_max_rt);
            if (p.rt_delta != 0) {
                ret += ", warped_rt: " +
                       std::to_string(p.local_max_rt + p.rt_delta);
                ret += " (" + std::to_string(p.rt_delta) + ")";
            }
            ret += ", local_max_height: " + std::to_string(p.local_max_height);
            ret += ", fitted_height: " + std::to_string(p.fitted_height);
            ret += ", fitted_mz: " + std::to_string(p.fitted_mz);
            ret += ", fitted_rt: " + std::to_string(p.fitted_rt);
            ret += ", fitted_sigma_mz: " + std::to_string(p.fitted_sigma_mz);
            ret += ", fitted_sigma_rt: " + std::to_string(p.fitted_sigma_rt);
            ret += ", fitted_volume: " + std::to_string(p.fitted_volume);
            ret += ">";
            return ret;
        });

    py::class_<Warp2D::TimeMap>(m, "TimeMap")
        .def_readonly("num_segments", &Warp2D::TimeMap::num_segments)
        .def_readonly("rt_start", &Warp2D::TimeMap::rt_start)
        .def_readonly("rt_end", &Warp2D::TimeMap::rt_end)
        .def_readonly("sample_rt_start", &Warp2D::TimeMap::sample_rt_start)
        .def_readonly("sample_rt_end", &Warp2D::TimeMap::sample_rt_end)
        .def("warp", &Warp2D::warp, py::arg("rt"))
        .def("__repr__", [](const Warp2D::TimeMap &m) {
            return "TimeMap <rt_min: " + std::to_string(m.rt_min) +
                   ", rt_max: " + std::to_string(m.rt_max) + ">";
        });

    py::class_<PythonAPI::SimilarityResults>(m, "Similarity")
        .def_readonly("self_a", &PythonAPI::SimilarityResults::self_a)
        .def_readonly("self_b", &PythonAPI::SimilarityResults::self_b)
        .def_readonly("overlap", &PythonAPI::SimilarityResults::overlap)
        .def_readonly("geometric_ratio",
                      &PythonAPI::SimilarityResults::geometric_ratio)
        .def_readonly("mean_ratio", &PythonAPI::SimilarityResults::mean_ratio)
        .def("__repr__", [](const PythonAPI::SimilarityResults &s) {
            return "Similarity: self_a: " + std::to_string(s.self_a) +
                   ", self_b: " + std::to_string(s.self_b) +
                   ", overlap: " + std::to_string(s.overlap) +
                   ", geometric_ratio: " + std::to_string(s.geometric_ratio) +
                   ", mean_ratio: " + std::to_string(s.mean_ratio);
        });

    py::class_<IdentData::SpectrumMatch>(m, "SpectrumMatch")
        .def_readonly("id", &IdentData::SpectrumMatch::id)
        .def_readonly("pass_threshold",
                      &IdentData::SpectrumMatch::pass_threshold)
        .def_readonly("match_id", &IdentData::SpectrumMatch::match_id)
        .def_readonly("charge_state", &IdentData::SpectrumMatch::charge_state)
        .def_readonly("theoretical_mz",
                      &IdentData::SpectrumMatch::theoretical_mz)
        .def_readonly("experimental_mz",
                      &IdentData::SpectrumMatch::experimental_mz)
        .def_readonly("retention_time",
                      &IdentData::SpectrumMatch::retention_time)
        .def_readonly("rank", &IdentData::SpectrumMatch::rank)
        .def("__repr__", [](const IdentData::SpectrumMatch &s) {
            return "SpectrumMatch <id: " + s.id +
                   ", pass_threshold: " + std::to_string(s.pass_threshold) +
                   ", match_id: " + s.match_id +
                   ", charge_state: " + std::to_string(s.charge_state) +
                   ", theoretical_mz: " + std::to_string(s.theoretical_mz) +
                   ", experimental_mz: " + std::to_string(s.experimental_mz) +
                   ", retention_time: " + std::to_string(s.retention_time) +
                   ", rank: " + std::to_string(s.rank) + ">";
        });

    py::class_<IdentData::DBSequence>(m, "DBSequence")
        .def_readonly("id", &IdentData::DBSequence::id)
        .def_readonly("accession", &IdentData::DBSequence::accession)
        .def_readonly("description", &IdentData::DBSequence::description)
        .def_readonly("db_reference", &IdentData::DBSequence::db_reference)
        .def("__repr__", [](const IdentData::DBSequence &s) {
            return "DBSequence <id: " + s.id + ", accession: " + s.accession +
                   ", description: " + s.description + ">";
        });

    py::class_<IdentData::PeptideModification>(m, "PeptideModification")
        .def_readonly("monoisotopic_mass_delta",
                      &IdentData::PeptideModification::monoisotopic_mass_delta)
        .def_readonly("average_mass_delta",
                      &IdentData::PeptideModification::average_mass_delta)
        .def_readonly("residues", &IdentData::PeptideModification::residues)
        .def_readonly("location", &IdentData::PeptideModification::location)
        .def_readonly("id", &IdentData::PeptideModification::id)
        .def("__repr__", [](const IdentData::PeptideModification &s) {
            return "PeptideModification <monoisotopic_mass_delta: " +
                   std::to_string(s.monoisotopic_mass_delta) +
                   ", average_mass_delta: " +
                   std::to_string(s.average_mass_delta) +
                   ", residues: " + s.residues +
                   ", location: " + std::to_string(s.location) + ">";
        });

    py::class_<IdentData::Peptide>(m, "Peptide")
        .def_readonly("id", &IdentData::Peptide::id)
        .def_readonly("sequence", &IdentData::Peptide::sequence)
        .def_readonly("modifications", &IdentData::Peptide::modifications)
        .def("__repr__", [](const IdentData::Peptide &s) {
            return "Peptide <id: " + s.id + ", sequence: " + s.sequence +
                   ", num_modifications: " +
                   std::to_string(s.modifications.size()) + ">";
        });

    py::class_<IdentData::PeptideEvidence>(m, "PeptideEvidence")
        .def_readonly("id", &IdentData::PeptideEvidence::id)
        .def_readonly("db_sequence_id",
                      &IdentData::PeptideEvidence::db_sequence_id)
        .def_readonly("peptide_id", &IdentData::PeptideEvidence::peptide_id)
        .def_readonly("decoy", &IdentData::PeptideEvidence::decoy)
        .def("__repr__", [](const IdentData::PeptideEvidence &s) {
            return "PeptideEvidence <id: " + s.id +
                   ", db_sequence_id: " + s.db_sequence_id +
                   ", peptide_id: " + s.peptide_id +
                   ", decoy: " + std::to_string(s.decoy) + ">";
        });

    py::class_<IdentData::IdentData>(m, "IdentData")
        .def_readonly("db_sequences", &IdentData::IdentData::db_sequences)
        .def_readonly("peptides", &IdentData::IdentData::peptides)
        .def_readonly("spectrum_matches",
                      &IdentData::IdentData::spectrum_matches)
        .def_readonly("peptide_evidence",
                      &IdentData::IdentData::peptide_evidence);

    py::class_<FeatureDetection::Feature>(m, "Feature")
        .def_readonly("id", &FeatureDetection::Feature::id)
        .def_readonly("score", &FeatureDetection::Feature::score)
        .def_readonly("average_rt", &FeatureDetection::Feature::average_rt)
        .def_readonly("average_rt_delta",
                      &FeatureDetection::Feature::average_rt_delta)
        .def_readonly("average_rt_sigma",
                      &FeatureDetection::Feature::average_rt_sigma)
        .def_readonly("average_mz", &FeatureDetection::Feature::average_mz)
        .def_readonly("average_mz_sigma",
                      &FeatureDetection::Feature::average_mz_sigma)
        .def_readonly("total_height", &FeatureDetection::Feature::total_height)
        .def_readonly("total_volume", &FeatureDetection::Feature::total_volume)
        .def_readonly("max_height", &FeatureDetection::Feature::max_height)
        .def_readonly("max_volume", &FeatureDetection::Feature::max_volume)
        .def_readonly("monoisotopic_mz",
                      &FeatureDetection::Feature::monoisotopic_mz)
        .def_readonly("monoisotopic_rt",
                      &FeatureDetection::Feature::monoisotopic_rt)
        .def_readonly("monoisotopic_height",
                      &FeatureDetection::Feature::monoisotopic_height)
        .def_readonly("monoisotopic_volume",
                      &FeatureDetection::Feature::monoisotopic_volume)
        .def_readonly("charge_state", &FeatureDetection::Feature::charge_state)
        .def_readonly("peak_ids", &FeatureDetection::Feature::peak_ids)
        .def("__repr__", [](const FeatureDetection::Feature &f) {
            std::string ret = "";
            ret += "Feature <id: " + std::to_string(f.id);
            ret += ", average_rt: " + std::to_string(f.average_rt);
            if (f.average_rt_delta != 0) {
                ret += ", average_warped_rt: " +
                       std::to_string(f.average_rt + f.average_rt_delta);
                ret += " (" + std::to_string(f.average_rt_delta) + ")";
            }
            ret += ", average_mz: " + std::to_string(f.average_mz);
            ret += ", total_height: " + std::to_string(f.total_height);
            ret += ", monoisotopic_mz: " + std::to_string(f.monoisotopic_mz);
            ret += ", monoisotopic_height: " +
                   std::to_string(f.monoisotopic_height);
            ret += ", charge_state: " + std::to_string(f.charge_state);
            ret += ", n_isotopes: " + std::to_string(f.peak_ids.size());
            ret += ">";
            return ret;
        });

    py::class_<MetaMatch::PeakId>(m, "PeakId")
        .def_readonly("file_id", &MetaMatch::PeakId::file_id)
        .def_readonly("peak_id", &MetaMatch::PeakId::peak_id)
        .def("__repr__", [](const MetaMatch::PeakId &c) {
            return "MetaClusterFileId <file_id: " + std::to_string(c.file_id) +
                   ", peak_id: " + std::to_string(c.peak_id) + ">";
        });

    py::class_<MetaMatch::FeatureId>(m, "FeatureId")
        .def_readonly("file_id", &MetaMatch::FeatureId::file_id)
        .def_readonly("feature_id", &MetaMatch::FeatureId::feature_id)
        .def("__repr__", [](const MetaMatch::FeatureId &c) {
            return "MetaClusterFileId <file_id: " + std::to_string(c.file_id) +
                   ", feature_id: " + std::to_string(c.feature_id) + ">";
        });

    py::class_<MetaMatch::PeakCluster>(m, "PeakCluster")
        .def_readonly("id", &MetaMatch::PeakCluster::id)
        .def_readonly("mz", &MetaMatch::PeakCluster::mz)
        .def_readonly("rt", &MetaMatch::PeakCluster::rt)
        .def_readonly("avg_height",
                      &MetaMatch::PeakCluster::avg_height)
        .def_readonly("avg_volume",
                      &MetaMatch::PeakCluster::avg_volume)
        .def_readonly("heights",
                      &MetaMatch::PeakCluster::heights)
        .def_readonly("volumes",
                      &MetaMatch::PeakCluster::volumes)
        .def_readonly("peak_ids", &MetaMatch::PeakCluster::peak_ids)
        .def("__repr__", [](const MetaMatch::PeakCluster &c) {
            return "MetaCluster <id: " + std::to_string(c.id) +
                   ", mz: " + std::to_string(c.mz) +
                   ", rt: " + std::to_string(c.rt) + ">";
        });

    py::class_<MetaMatch::FeatureCluster>(m, "FeatureCluster")
        .def_readonly("id", &MetaMatch::FeatureCluster::id)
        .def_readonly("mz", &MetaMatch::FeatureCluster::mz)
        .def_readonly("rt", &MetaMatch::FeatureCluster::rt)
        .def_readonly("avg_total_height",
                      &MetaMatch::FeatureCluster::avg_total_height)
        .def_readonly("avg_monoisotopic_height",
                      &MetaMatch::FeatureCluster::avg_monoisotopic_height)
        .def_readonly("avg_max_height",
                      &MetaMatch::FeatureCluster::avg_max_height)
        .def_readonly("avg_total_volume",
                      &MetaMatch::FeatureCluster::avg_total_volume)
        .def_readonly("avg_monoisotopic_volume",
                      &MetaMatch::FeatureCluster::avg_monoisotopic_volume)
        .def_readonly("avg_max_volume",
                      &MetaMatch::FeatureCluster::avg_max_volume)
        .def_readonly("charge_state", &MetaMatch::FeatureCluster::charge_state)
        .def_readonly("total_heights",
                      &MetaMatch::FeatureCluster::total_heights)
        .def_readonly("monoisotopic_heights",
                      &MetaMatch::FeatureCluster::monoisotopic_heights)
        .def_readonly("max_heights", &MetaMatch::FeatureCluster::max_heights)
        .def_readonly("total_volumes",
                      &MetaMatch::FeatureCluster::total_volumes)
        .def_readonly("monoisotopic_volumes",
                      &MetaMatch::FeatureCluster::monoisotopic_volumes)
        .def_readonly("max_volumes", &MetaMatch::FeatureCluster::max_volumes)
        .def_readonly("feature_ids", &MetaMatch::FeatureCluster::feature_ids)
        .def("__repr__", [](const MetaMatch::FeatureCluster &c) {
            return "MetaCluster <id: " + std::to_string(c.id) +
                   ", mz: " + std::to_string(c.mz) +
                   ", rt: " + std::to_string(c.rt) + ">";
        });

    py::class_<Link::LinkedMsms>(m, "LinkedMsms")
        .def_readonly("entity_id", &Link::LinkedMsms::entity_id)
        .def_readonly("msms_id", &Link::LinkedMsms::msms_id)
        .def_readonly("distance", &Link::LinkedMsms::distance)
        .def("__repr__", [](const Link::LinkedMsms &p) {
            return "LinkedMsms <entity_id: " + std::to_string(p.entity_id) +
                   ", scan_index: " + std::to_string(p.scan_index) +
                   ", msms_id: " + std::to_string(p.msms_id) +
                   ", distance: " + std::to_string(p.distance) + ">";
        });

    py::class_<Link::LinkedPsm>(m, "LinkedPsm")
        .def_readonly("peak_id", &Link::LinkedPsm::peak_id)
        .def_readonly("psm_index", &Link::LinkedPsm::psm_index)
        .def_readonly("distance", &Link::LinkedPsm::distance)
        .def("__repr__", [](const Link::LinkedPsm &p) {
            return "LinkedPsm <peak_id: " + std::to_string(p.peak_id) +
                   ", psm_index: " + std::to_string(p.psm_index) +
                   ", distance: " + std::to_string(p.distance) + ">";
        });


    py::class_<ProteinInference::InferredProtein>(m, "InferredProtein")
        .def_readonly("protein_id",
                      &ProteinInference::InferredProtein::protein_id)
        .def_readonly("psm_id", &ProteinInference::InferredProtein::psm_id)
        .def("__repr__", [](const ProteinInference::InferredProtein &p) {
            return "InferredProtein <protein_id: " + p.protein_id +
                   ", psm_id: " + p.psm_id + ">";
        });

    // Functions.
    m.def("read_mzxml", &PythonAPI::read_mzxml,
          "Read raw data from the given mzXML file ", py::arg("file_name"),
          py::arg("min_mz") = -1.0, py::arg("max_mz") = -1.0,
          py::arg("min_rt") = -1.0, py::arg("max_rt") = -1.0,
          py::arg("instrument_type") = "", py::arg("resolution_ms1"),
          py::arg("resolution_msn"), py::arg("reference_mz"),
          py::arg("fwhm_rt"), py::arg("polarity") = "", py::arg("ms_level") = 1)
        .def("read_mzml", &PythonAPI::read_mzml,
             "Read raw data from the given mzXML file ", py::arg("file_name"),
             py::arg("min_mz") = -1.0, py::arg("max_mz") = -1.0,
             py::arg("min_rt") = -1.0, py::arg("max_rt") = -1.0,
             py::arg("instrument_type") = "", py::arg("resolution_ms1"),
             py::arg("resolution_msn"), py::arg("reference_mz"),
             py::arg("fwhm_rt"), py::arg("polarity") = "",
             py::arg("ms_level") = 1)
        .def("theoretical_fwhm", &RawData::theoretical_fwhm,
             "Calculate the theoretical width of the peak at the given m/z for "
             "the given raw file",
             py::arg("raw_data"), py::arg("mz"))
        .def("resample", &PythonAPI::resample,
             "Resample the raw data into a smoothed warped grid",
             py::arg("raw_data"), py::arg("num_mz") = 10,
             py::arg("num_rt") = 10, py::arg("smoothing_coef_mz") = 0.5,
             py::arg("smoothing_coef_rt") = 0.5)
        .def("find_peaks", &Centroid::find_peaks_parallel,
             "Find all peaks in the given grid", py::arg("raw_data"),
             py::arg("grid"), py::arg("max_peaks") = 0,
             py::arg("max_threads") = std::thread::hardware_concurrency())
        .def("calculate_time_map", &PythonAPI::calculate_time_map,
             "Calculate a warping time_map to maximize the similarity of "
             "ref_peaks and source_peaks",
             py::arg("ref_peaks"), py::arg("source_peaks"), py::arg("slack"),
             py::arg("window_size"), py::arg("num_points"),
             py::arg("rt_expand_factor"), py::arg("peaks_per_window"))
        .def("warp_peaks", &Warp2D::warp_peaks,
             "Warp the peak list using the given time map", py::arg("peaks"),
             py::arg("time_map"))
        .def("find_similarity", &PythonAPI::find_similarity,
             "Find the similarity between two peak lists",
             py::arg("peak_list_a"), py::arg("peak_list_b"), py::arg("n_peaks"))
        .def("write_peaks", &PythonAPI::write_peaks,
             "Write the peaks to disk in a binary format", py::arg("peaks"),
             py::arg("file_name"))
        .def("write_linked_msms", &PythonAPI::write_linked_msms,
             "Write the linked_msms to disk in a binary format",
             py::arg("linked_msms"), py::arg("file_name"))
        .def("write_linked_psm", &PythonAPI::write_linked_psm,
             "Write the linked_psm to disk in a binary format",
             py::arg("linked_psm"), py::arg("file_name"))
        .def("read_peaks", &PythonAPI::read_peaks,
             "Read the peaks from the binary peaks file", py::arg("file_name"))
        .def("read_raw_data", &PythonAPI::read_raw_data,
             "Read the raw_data from the binary raw_data file",
             py::arg("file_name"))
        .def("read_grid", &PythonAPI::read_grid,
             "Read the grid from the binary grid file", py::arg("file_name"))
        .def("read_linked_psm", &PythonAPI::read_linked_psm,
             "Read the linked_psm from the binary linked_psm file",
             py::arg("file_name"))
        .def("read_linked_msms", &PythonAPI::read_linked_msms,
             "Read the linked_msms from the binary linked_msms file",
             py::arg("file_name"))
        .def("read_time_map", &PythonAPI::read_time_map,
             "Read the time_map from the binary time_map file",
             py::arg("file_name"))
        .def("write_time_map", &PythonAPI::write_time_map,
             "Write the time_map to disk in a binary format",
             py::arg("time_map"), py::arg("file_name"))
        .def("read_mzidentml", &PythonAPI::read_mzidentml,
             "Read identification data from the given mzIdentML file ",
             py::arg("file_name"), py::arg("ignore_decoy") = true,
             py::arg("require_threshold") = true,
             py::arg("max_rank_only") = true,
             py::arg("min_mz") = -1.0, py::arg("max_mz") = -1.0,
             py::arg("min_rt") = -1.0, py::arg("max_rt") = -1.0)
        .def("read_ident_data", &PythonAPI::read_ident_data,
             "Read the ident_data from the binary ident_data file",
             py::arg("file_name"))
        .def("write_ident_data", &PythonAPI::write_ident_data,
             "Write the ident_data to disk in a binary format",
             py::arg("ident_data"), py::arg("file_name"))
        .def(
            "read_inferred_proteins", &PythonAPI::read_inferred_proteins,
            "Read the inferred_proteins from the binary inferred_proteins file",
            py::arg("file_name"))
        .def("write_inferred_proteins", &PythonAPI::write_inferred_proteins,
             "Write the inferred_proteins to disk in a binary format",
             py::arg("inferred_proteins"), py::arg("file_name"))
        .def("read_features", &PythonAPI::read_features,
             "Read the feature from the binary feature file",
             py::arg("file_name"))
        .def("write_features", &PythonAPI::write_features,
             "Write the feature to disk in a binary format", py::arg("feature"),
             py::arg("file_name"))
        .def("read_peak_clusters", &PythonAPI::read_peak_clusters,
             "Read the peak_clusters from the binary peak_clusters file",
             py::arg("file_name"))
        .def("write_peak_clusters", &PythonAPI::write_peak_clusters,
             "Write the peak_clusters to disk in a binary format",
             py::arg("peak_clusters"), py::arg("file_name"))
        .def("read_feature_clusters", &PythonAPI::read_feature_clusters,
             "Read the feature_clusters from the binary feature_clusters file",
             py::arg("file_name"))
        .def("write_feature_clusters", &PythonAPI::write_feature_clusters,
             "Write the feature_clusters to disk in a binary format",
             py::arg("feature_clusters"), py::arg("file_name"))
        .def("find_feature_clusters", &PythonAPI::find_feature_clusters,
             "Perform metamatch for feature matching", py::arg("group_ids"),
             py::arg("features"), py::arg("keep_perc"),
             py::arg("intensity_threshold") = 0.5, py::arg("n_sig_mz") = 1.5,
             py::arg("n_sig_rt") = 1.5)
        .def("find_peak_clusters", &PythonAPI::find_peak_clusters,
             "Perform metamatch for peak matching", py::arg("group_ids"),
             py::arg("features"), py::arg("keep_perc"),
             py::arg("intensity_threshold") = 0.5, py::arg("n_sig_mz") = 1.5,
             py::arg("n_sig_rt") = 1.5)
        .def("link_peaks", &Link::link_peaks, "Link msms events to peak ids",
             py::arg("peaks"), py::arg("raw_data"), py::arg("n_sig_mz") = 3,
             py::arg("n_sig_rt") = 3)
        .def("link_idents", &Link::link_idents,
             "Link msms events to spectrum identifications",
             py::arg("ident_data"), py::arg("raw_data"),
             py::arg("n_sig_mz") = 3, py::arg("n_sig_rt") = 3)
        .def("link_psm", &Link::link_psm,
             "Link spectrum identifications with peaks",
             py::arg("ident_data"), py::arg("peaks"), py::arg("raw_data"),
             py::arg("n_sig_mz") = 3, py::arg("n_sig_rt") = 3)
        .def("xic", &PythonAPI::xic, py::arg("raw_data"), py::arg("min_mz"),
             py::arg("max_mz"), py::arg("min_rt"), py::arg("max_rt"),
             py::arg("method") = "sum")
        .def("perform_protein_inference", &ProteinInference::razor,
             py::arg("ident_data"))
        .def("detect_features", &FeatureDetection::detect_features,
             "Link peaks as features", py::arg("peaks"),
             py::arg("charge_states"));
}
