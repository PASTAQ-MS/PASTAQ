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
#include "raw_data/raw_data.hpp"
#include "raw_data/raw_data_serialize.hpp"
#include "raw_data/xml_reader.hpp"
#include "utils/search.hpp"
#include "utils/serialization.hpp"
#include "warp2d/warp2d.hpp"

namespace py = pybind11;

namespace PythonAPI {
RawData::RawData read_mzxml(std::string &input_file, double min_mz,
                            double max_mz, double min_rt, double max_rt,
                            std::string instrument_type_str,
                            double resolution_ms1, double resolution_msn,
                            double reference_mz, double fwhm_rt,
                            std::string polarity_str, size_t ms_level) {
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
        std::ostringstream error_stream;
        error_stream << "the given polarity is not supported. choose "
                        "between '+', '-', 'both' (default)";
        throw std::invalid_argument(error_stream.str());
    }

    // Sanity check the min/max rt/mz.
    if (min_rt >= max_rt) {
        std::ostringstream error_stream;
        error_stream << "error: min_rt >= max_rt (min_rt: " << min_rt
                     << ", max_rt: " << max_rt << ")";
        throw std::invalid_argument(error_stream.str());
    }
    if (min_mz >= max_mz) {
        std::ostringstream error_stream;
        error_stream << "error: min_mz >= max_mz (min_mz: " << min_mz
                     << ", max_mz: " << max_mz << ")";
        throw std::invalid_argument(error_stream.str());
    }

    // Open file stream.
    std::ifstream stream;
    stream.open(input_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }

    auto raw_data = XmlReader::read_mzxml(
        stream, min_mz, max_mz, min_rt, max_rt, instrument_type, resolution_ms1,
        resolution_msn, reference_mz, polarity, ms_level);
    if (!raw_data) {
        std::ostringstream error_stream;
        error_stream << "error: an error occurred when reading the file"
                     << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    raw_data.value().fwhm_rt = fwhm_rt;

    return raw_data.value();
}

Xic::Xic xic(const RawData::RawData &raw_data, double min_mz, double max_mz,
             double min_rt, double max_rt, std::string method_str) {
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
        std::ostringstream error_stream;
        error_stream << "the given xic method is not supported";
        throw std::invalid_argument(error_stream.str());
    }
    return RawData::xic(raw_data, min_mz, max_mz, min_rt, max_rt, method);
}

Grid::Grid resample(const RawData::RawData &raw_data, uint64_t num_samples_mz,
                    uint64_t num_samples_rt, double smoothing_coef_mz,
                    double smoothing_coef_rt) {
    auto params = Grid::ResampleParams{};
    params.num_samples_mz = num_samples_mz;
    params.num_samples_rt = num_samples_rt;
    params.smoothing_coef_mz = smoothing_coef_mz;
    params.smoothing_coef_rt = smoothing_coef_rt;
    return Grid::resample(raw_data, params);
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

std::vector<std::vector<Centroid::Peak>> warp_peaks(
    const std::vector<std::vector<Centroid::Peak>> &all_peaks,
    size_t reference_index, int64_t slack, int64_t window_size,
    int64_t num_points, double rt_expand_factor, int64_t peaks_per_window) {
    // TODO(alex): Validate the parameters and throw an error if
    // appropriate.
    Warp2D::Parameters parameters = {slack, window_size, num_points,
                                     peaks_per_window, rt_expand_factor};
    auto reference_peaks = all_peaks[reference_index];

    std::vector<std::vector<Centroid::Peak>> all_warped_peaks;
    for (size_t i = 0; i < all_peaks.size(); ++i) {
        if (i == reference_index) {
            all_warped_peaks.push_back(all_peaks[i]);
            continue;
        }
        auto peaks = all_peaks[i];
        std::vector<Centroid::Peak> warped_peaks;
        warped_peaks =
            Warp2D::warp_peaks_parallel(reference_peaks, peaks, parameters,
                                        std::thread::hardware_concurrency());
        all_warped_peaks.push_back(warped_peaks);
    }
    return all_warped_peaks;
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
    auto sort_peaks = [](const Centroid::Peak &p1,
                         const Centroid::Peak &p2) -> bool {
        return (p1.local_max_height >= p2.local_max_height);
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
    // Overlap / (GeometricMean(self_a, self_b))
    results.geometric_ratio =
        results.overlap / std::sqrt(results.self_a * results.self_b);
    // Harmonic mean of the ratios between
    // self_similarity/overlap_similarity
    results.mean_ratio =
        2 * results.overlap / (results.self_a + results.self_b);
    return results;
}

void write_raw_data(const RawData::RawData &raw_data,
                    std::string &output_file) {
    // Open file stream.
    std::ofstream stream;
    stream.open(output_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open output file" << output_file;
        throw std::invalid_argument(error_stream.str());
    }

    if (!RawData::Serialize::write_raw_data(stream, raw_data)) {
        std::ostringstream error_stream;
        error_stream
            << "error: couldn't write the raw_data into the output file"
            << output_file;
        throw std::invalid_argument(error_stream.str());
    }
}

RawData::RawData read_raw_data(std::string &input_file) {
    // Open file stream.
    std::ifstream stream;
    stream.open(input_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }

    RawData::RawData raw_data;
    if (!RawData::Serialize::read_raw_data(stream, &raw_data)) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't write the raw_data into the input file"
                     << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    return raw_data;
}

void write_grid(const Grid::Grid &grid, std::string &output_file) {
    // Open file stream.
    std::ofstream stream;
    stream.open(output_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open output file" << output_file;
        throw std::invalid_argument(error_stream.str());
    }

    if (!Grid::Serialize::write_grid(stream, grid)) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't write the grid into the output file"
                     << output_file;
        throw std::invalid_argument(error_stream.str());
    }
}

Grid::Grid read_grid(std::string &input_file) {
    // Open file stream.
    std::ifstream stream;
    stream.open(input_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }

    Grid::Grid grid;
    if (!Grid::Serialize::read_grid(stream, &grid)) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't write the grid into the input file"
                     << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    return grid;
}

void write_peaks(const std::vector<Centroid::Peak> &peaks,
                 std::string &output_file) {
    // Open file stream.
    std::ofstream stream;
    stream.open(output_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open output file" << output_file;
        throw std::invalid_argument(error_stream.str());
    }

    if (!Centroid::Serialize::write_peaks(stream, peaks)) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't write the peaks into the output file"
                     << output_file;
        throw std::invalid_argument(error_stream.str());
    }
}

std::vector<Centroid::Peak> read_peaks(std::string &input_file) {
    // Open file stream.
    std::ifstream stream;
    stream.open(input_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }

    std::vector<Centroid::Peak> peaks;
    if (!Centroid::Serialize::read_peaks(stream, &peaks)) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't write the peaks into the input file"
                     << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    return peaks;
}

void write_feature_clusters(
    const std::vector<MetaMatch::FeatureCluster> &feature_clusters,
    std::string &output_file) {
    // Open file stream.
    std::ofstream stream;
    stream.open(output_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open output file" << output_file;
        throw std::invalid_argument(error_stream.str());
    }

    if (!MetaMatch::Serialize::write_feature_clusters(stream,
                                                      feature_clusters)) {
        std::ostringstream error_stream;
        error_stream
            << "error: couldn't write the feature_clusters into the output file"
            << output_file;
        throw std::invalid_argument(error_stream.str());
    }
}

std::vector<MetaMatch::FeatureCluster> read_feature_clusters(
    std::string &input_file) {
    // Open file stream.
    std::ifstream stream;
    stream.open(input_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }

    std::vector<MetaMatch::FeatureCluster> feature_clusters;
    if (!MetaMatch::Serialize::read_feature_clusters(stream,
                                                     &feature_clusters)) {
        std::ostringstream error_stream;
        error_stream
            << "error: couldn't write the feature_clusters into the input file"
            << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    return feature_clusters;
}

void write_features(const std::vector<FeatureDetection::Feature> &features,
                    std::string &output_file) {
    // Open file stream.
    std::ofstream stream;
    stream.open(output_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open output file" << output_file;
        throw std::invalid_argument(error_stream.str());
    }

    if (!FeatureDetection::Serialize::write_features(stream, features)) {
        std::ostringstream error_stream;
        error_stream
            << "error: couldn't write the features into the output file"
            << output_file;
        throw std::invalid_argument(error_stream.str());
    }
}

std::vector<FeatureDetection::Feature> read_features(std::string &input_file) {
    // Open file stream.
    std::ifstream stream;
    stream.open(input_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }

    std::vector<FeatureDetection::Feature> features;
    if (!FeatureDetection::Serialize::read_features(stream, &features)) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't write the features into the input file"
                     << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    return features;
}

void write_ident_data(const IdentData::IdentData &ident_data,
                      std::string &output_file) {
    // Open file stream.
    std::ofstream stream;
    stream.open(output_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open output file" << output_file;
        throw std::invalid_argument(error_stream.str());
    }

    if (!IdentData::Serialize::write_ident_data(stream, ident_data)) {
        std::ostringstream error_stream;
        error_stream
            << "error: couldn't write the ident_data into the output file"
            << output_file;
        throw std::invalid_argument(error_stream.str());
    }
}

IdentData::IdentData read_ident_data(std::string &input_file) {
    // Open file stream.
    std::ifstream stream;
    stream.open(input_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }

    IdentData::IdentData ident_data;
    if (!IdentData::Serialize::read_ident_data(stream, &ident_data)) {
        std::ostringstream error_stream;
        error_stream
            << "error: couldn't write the ident_data into the input file"
            << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    return ident_data;
}

void write_metamatch_clusters(
    const std::vector<MetaMatch::Cluster> &metamatch_clusters,
    std::string &output_file) {
    // Open file stream.
    std::ofstream stream;
    stream.open(output_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open output file" << output_file;
        throw std::invalid_argument(error_stream.str());
    }

    if (!MetaMatch::Serialize::write_clusters(stream, metamatch_clusters)) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't write the metamatch_clusters into the "
                        "output file"
                     << output_file;
        throw std::invalid_argument(error_stream.str());
    }
}

std::vector<MetaMatch::Cluster> read_metamatch_clusters(
    std::string &input_file) {
    // Open file stream.
    std::ifstream stream;
    stream.open(input_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }

    std::vector<MetaMatch::Cluster> metamatch_clusters;
    if (!MetaMatch::Serialize::read_clusters(stream, &metamatch_clusters)) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't write the metamatch_clusters into the "
                        "input file"
                     << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    return metamatch_clusters;
}

void write_metamatch_peaks(const std::vector<MetaMatch::Peak> &metamatch_peaks,
                           std::string &output_file) {
    // Open file stream.
    std::ofstream stream;
    stream.open(output_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open output file" << output_file;
        throw std::invalid_argument(error_stream.str());
    }

    if (!MetaMatch::Serialize::write_peaks(stream, metamatch_peaks)) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't write the metamatch_peaks into the "
                        "output file"
                     << output_file;
        throw std::invalid_argument(error_stream.str());
    }
}

std::vector<MetaMatch::Peak> read_metamatch_peaks(std::string &input_file) {
    // Open file stream.
    std::ifstream stream;
    stream.open(input_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }

    std::vector<MetaMatch::Peak> metamatch_peaks;
    if (!MetaMatch::Serialize::read_peaks(stream, &metamatch_peaks)) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't write the metamatch_peaks into the "
                        "input file"
                     << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    return metamatch_peaks;
}

void write_linked_msms(const std::vector<Link::LinkedMsms> &linked_msms,
                       std::string &output_file) {
    // Open file stream.
    std::ofstream stream;
    stream.open(output_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open output file" << output_file;
        throw std::invalid_argument(error_stream.str());
    }

    if (!Link::Serialize::write_linked_msms_table(stream, linked_msms)) {
        std::ostringstream error_stream;
        error_stream
            << "error: couldn't write the linked_msms into the output file"
            << output_file;
        throw std::invalid_argument(error_stream.str());
    }
}

std::vector<Link::LinkedMsms> read_linked_msms(std::string &input_file) {
    // Open file stream.
    std::ifstream stream;
    stream.open(input_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }

    std::vector<Link::LinkedMsms> linked_msms;
    if (!Link::Serialize::read_linked_msms_table(stream, &linked_msms)) {
        std::ostringstream error_stream;
        error_stream
            << "error: couldn't write the linked_msms into the input file"
            << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    return linked_msms;
}

IdentData::IdentData read_mzidentml(std::string &input_file) {
    // Open file stream.
    std::ifstream stream;
    stream.open(input_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    return XmlReader::read_mzidentml(stream);
}

struct MetaMatchResults {
    std::vector<MetaMatch::Cluster> clusters;
    std::vector<MetaMatch::Peak> orphans;
};

MetaMatchResults perform_metamatch(
    // NOTE: [(class_0, peaks_0),...(class_i, peaks_i)]
    std::vector<std::tuple<uint32_t, std::vector<Centroid::Peak>>> input,
    double radius_mz, double radius_rt, double fraction) {
    MetaMatchResults results;

    // Create the ClassMaps.
    std::vector<MetaMatch::ClassMap> class_maps;
    std::vector<MetaMatch::Peak> metapeaks;
    uint32_t file_id = 0;
    for (const auto &[class_id, peaks] : input) {
        bool found = false;
        for (auto &class_map : class_maps) {
            if (class_map.id == class_id) {
                found = true;
                class_map.n_files++;
                class_map.required_hits =
                    std::ceil(class_map.n_files * fraction);
                break;
            }
        }
        if (!found) {
            class_maps.push_back({class_id, 1, 0});
        }
        for (const auto &peak : peaks) {
            metapeaks.push_back({peak, file_id, class_id, -1, peak.local_max_mz,
                                 peak.local_max_rt});
        }
        ++file_id;
    }
    MetaMatch::Parameters parameters = {};
    parameters.radius_mz = radius_mz;
    parameters.radius_rt = radius_rt;
    parameters.class_maps = class_maps;

    MetaMatch::find_clusters(metapeaks, parameters);
    results.orphans = MetaMatch::extract_orphans(metapeaks);
    results.clusters = MetaMatch::reduce_cluster(metapeaks, file_id);

    return results;
}

std::vector<MetaMatch::FeatureCluster> find_feature_clusters(
    std::vector<uint64_t> group_ids,
    std::vector<std::vector<Centroid::Peak>> peaks,
    std::vector<std::vector<FeatureDetection::Feature>> features) {
    if (group_ids.size() != peaks.size() ||
        group_ids.size() != features.size()) {
        std::ostringstream error_stream;
        error_stream
            << "error: groups, peaks and features have different lengths";
        throw std::invalid_argument(error_stream.str());
    }
    // Create input set.
    std::vector<MetaMatch::InputSetFeatures> input_sets;
    for (size_t i = 0; i < group_ids.size(); ++i) {
        MetaMatch::InputSetFeatures input_set = {group_ids[i], peaks[i],
                                                 features[i]};
        input_sets.push_back(input_set);
    }
    return MetaMatch::find_feature_clusters(input_sets);
}

void debug(IdentData::IdentData &ident_data) {
    auto inference_graph = ProteinInference::create_graph(ident_data);
    ProteinInference::razor(inference_graph);
}

}  // namespace PythonAPI

PYBIND11_MODULE(tapp, m) {
    // Documentation.
    m.doc() = "tapp documentation";

    // Structs.
    py::class_<RawData::Scan>(m, "Scan")
        .def_readonly("scan_number", &RawData::Scan::scan_number)
        .def_readonly("ms_level", &RawData::Scan::ms_level)
        .def_readonly("num_points", &RawData::Scan::num_points)
        .def_readonly("retention_time", &RawData::Scan::retention_time)
        .def_readonly("mz", &RawData::Scan::mz)
        .def_readonly("intensity", &RawData::Scan::intensity)
        .def_readonly("polarity", &RawData::Scan::polarity)
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
            ret += ", raw_roi_sigma_mz: " + std::to_string(p.raw_roi_sigma_mz);
            ret += ", raw_roi_sigma_rt: " + std::to_string(p.raw_roi_sigma_rt);
            ret +=
                ", raw_roi_num_points: " + std::to_string(p.raw_roi_num_points);
            ret +=
                ", raw_roi_num_scans: " + std::to_string(p.raw_roi_num_scans);
            ret += ">";
            return ret;
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

    py::class_<IdentData::SpectrumId>(m, "SpectrumId")
        .def_readonly("id", &IdentData::SpectrumId::id)
        .def_readonly("pass_threshold", &IdentData::SpectrumId::pass_threshold)
        .def_readonly("modifications", &IdentData::SpectrumId::modifications)
        .def_readonly("sequence", &IdentData::SpectrumId::sequence)
        .def_readonly("peptide_id", &IdentData::SpectrumId::peptide_id)
        .def_readonly("charge_state", &IdentData::SpectrumId::charge_state)
        .def_readonly("theoretical_mz", &IdentData::SpectrumId::theoretical_mz)
        .def_readonly("experimental_mz",
                      &IdentData::SpectrumId::experimental_mz)
        .def_readonly("retention_time", &IdentData::SpectrumId::retention_time)
        .def_readonly("rank", &IdentData::SpectrumId::rank)
        .def("__repr__", [](const IdentData::SpectrumId &s) {
            return "SpectrumId <id: " + s.id +
                   ", pass_threshold: " + std::to_string(s.pass_threshold) +
                   ", modifications: " + std::to_string(s.modifications) +
                   ", sequence: " + s.sequence +
                   ", peptide_id: " + s.peptide_id +
                   ", charge_state: " + std::to_string(s.charge_state) +
                   ", theoretical_mz: " + std::to_string(s.theoretical_mz) +
                   ", experimental_mz: " + std::to_string(s.experimental_mz) +
                   ", retention_time: " + std::to_string(s.retention_time) +
                   ", rank: " + std::to_string(s.rank) + ">";
        });

    py::class_<IdentData::DBSequence>(m, "DBSequence")
        .def_readonly("id", &IdentData::DBSequence::id)
        .def_readonly("value", &IdentData::DBSequence::value)
        .def("__repr__", [](const IdentData::DBSequence &s) {
            return "DBSequence <id: " + s.id + ", value: " + s.value + ">";
        });

    py::class_<IdentData::PeptideModification>(m, "PeptideModification")
        .def_readonly("monoisotopic_mass_delta",
                      &IdentData::PeptideModification::monoisotopic_mass_delta)
        .def_readonly("average_mass_delta",
                      &IdentData::PeptideModification::average_mass_delta)
        .def_readonly("residues", &IdentData::PeptideModification::residues)
        .def_readonly("location", &IdentData::PeptideModification::location)
        .def("__repr__", [](const IdentData::PeptideModification &s) {
            return "PeptideModification <monoisotopic_mass_delta: " +
                   std::to_string(s.monoisotopic_mass_delta) +
                   ", average_mass_delta: " +
                   std::to_string(s.average_mass_delta) +
                   ", residues: " + s.residues +
                   ", location: " + std::to_string(s.location) +
                   ", num_cv_params: " + std::to_string(s.cv_params.size()) +
                   ">";
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

    py::class_<IdentData::ProteinHypothesis>(m, "ProteinHypothesis")
        .def_readonly("db_sequence_id",
                      &IdentData::ProteinHypothesis::db_sequence_id)
        .def_readonly("pass_threshold",
                      &IdentData::ProteinHypothesis::pass_threshold)
        .def_readonly("spectrum_ids",
                      &IdentData::ProteinHypothesis::spectrum_ids)
        .def("__repr__", [](const IdentData::ProteinHypothesis &s) {
            return "ProteinHypothesis <db_sequence_id: " + s.db_sequence_id +
                   ", pass_threshold: " + std::to_string(s.pass_threshold) +
                   ", num_spectrum_ids: " +
                   std::to_string(s.spectrum_ids.size()) + ">";
        });

    py::class_<IdentData::IdentData>(m, "IdentData")
        .def_readonly("db_sequences", &IdentData::IdentData::db_sequences)
        .def_readonly("peptides", &IdentData::IdentData::peptides)
        .def_readonly("spectrum_ids", &IdentData::IdentData::spectrum_ids)
        .def_readonly("protein_hypotheses",
                      &IdentData::IdentData::protein_hypotheses);

    py::class_<FeatureDetection::Feature>(m, "Feature")
        .def_readonly("id", &FeatureDetection::Feature::id)
        .def_readonly("msms_id", &FeatureDetection::Feature::msms_id)
        .def_readonly("average_rt", &FeatureDetection::Feature::average_rt)
        .def_readonly("average_rt_delta",
                      &FeatureDetection::Feature::average_rt_delta)
        .def_readonly("average_rt_sigma",
                      &FeatureDetection::Feature::average_rt_sigma)
        .def_readonly("average_mz", &FeatureDetection::Feature::average_mz)
        .def_readonly("average_mz_sigma",
                      &FeatureDetection::Feature::average_mz_sigma)
        .def_readonly("total_height", &FeatureDetection::Feature::total_height)
        .def_readonly("monoisotopic_mz",
                      &FeatureDetection::Feature::monoisotopic_mz)
        .def_readonly("monoisotopic_height",
                      &FeatureDetection::Feature::monoisotopic_height)
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
            ret += ", n_isotopes: " + std::to_string(f.peak_ids.size());
            ret += ">";
            return ret;
        });

    py::class_<PythonAPI::MetaMatchResults>(m, "MetaMatchResults")
        .def_readonly("clusters", &PythonAPI::MetaMatchResults::clusters)
        .def_readonly("orphans", &PythonAPI::MetaMatchResults::orphans);

    py::class_<MetaMatch::Cluster>(m, "MetaMatchCluster")
        .def_readonly("id", &MetaMatch::Cluster::id)
        .def_readonly("mz", &MetaMatch::Cluster::mz)
        .def_readonly("rt", &MetaMatch::Cluster::rt)
        .def_readonly("file_heights", &MetaMatch::Cluster::file_heights)
        .def_readonly("avg_height", &MetaMatch::Cluster::avg_height)
        .def("__repr__", [](const MetaMatch::Cluster &c) {
            return "MetaCluster <id: " + std::to_string(c.id) +
                   ", mz: " + std::to_string(c.mz) +
                   ", rt: " + std::to_string(c.rt) +
                   ", avg_height: " + std::to_string(c.avg_height) + ">";
        });

    py::class_<MetaMatch::FeatureCluster>(m, "FeatureCluster")
        .def_readonly("id", &MetaMatch::FeatureCluster::id)
        .def_readonly("mz", &MetaMatch::FeatureCluster::mz)
        .def_readonly("rt", &MetaMatch::FeatureCluster::rt)
        .def_readonly("avg_height", &MetaMatch::FeatureCluster::avg_height)
        .def_readonly("file_heights", &MetaMatch::FeatureCluster::file_heights)
        .def("__repr__", [](const MetaMatch::FeatureCluster &c) {
            return "MetaCluster <id: " + std::to_string(c.id) +
                   ", mz: " + std::to_string(c.mz) +
                   ", rt: " + std::to_string(c.rt) +
                   ", avg_height: " + std::to_string(c.avg_height) + ">";
        });

    py::class_<MetaMatch::Peak>(m, "MetaMatchPeak")
        .def_readonly("file_id", &MetaMatch::Peak::file_id)
        .def_readonly("class_id", &MetaMatch::Peak::class_id)
        .def_readonly("cluster_id", &MetaMatch::Peak::cluster_id)
        .def_readonly("cluster_mz", &MetaMatch::Peak::cluster_mz)
        .def_readonly("cluster_rt", &MetaMatch::Peak::cluster_rt)
        .def_readonly("height", &MetaMatch::Peak::local_max_height)
        .def_readonly("local_max_mz", &MetaMatch::Peak::local_max_mz)
        .def_readonly("local_max_rt", &MetaMatch::Peak::local_max_rt)
        .def("__repr__", [](const MetaMatch::Peak &p) {
            return "MetaPeak <peak_id: " + std::to_string(p.id) +
                   ", file_id: " + std::to_string(p.file_id) +
                   ", class_id: " + std::to_string(p.class_id) + ">";
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

    // Functions.
    m.def("read_mzxml", &PythonAPI::read_mzxml,
          "Read raw data from the given mzXML file ", py::arg("file_name"),
          py::arg("min_mz") = -1.0, py::arg("max_mz") = -1.0,
          py::arg("min_rt") = -1.0, py::arg("max_rt") = -1.0,
          py::arg("instrument_type") = "", py::arg("resolution_ms1"),
          py::arg("resolution_msn"), py::arg("reference_mz"),
          py::arg("fwhm_rt"), py::arg("polarity") = "", py::arg("ms_level") = 1)
        .def("theoretical_fwhm", &RawData::theoretical_fwhm,
             "Calculate the theoretical width of the peak at the given m/z for "
             "the given raw file",
             py::arg("raw_data"), py::arg("mz"))
        .def("resample", &PythonAPI::resample,
             "Resample the raw data into a smoothed warped grid",
             py::arg("raw_data"), py::arg("num_mz") = 10,
             py::arg("num_rt") = 10, py::arg("smoothing_coef_mz") = 0.5,
             py::arg("smoothing_coef_rt") = 0.5)
        .def("find_raw_points", &RawData::find_raw_points,
             "Save the fitted peaks as a bpks file", py::arg("raw_data"),
             py::arg("min_mz"), py::arg("max_mz"), py::arg("min_rt"),
             py::arg("max_rt"))
        .def("find_peaks", &Centroid::find_peaks_parallel,
             "Find all peaks in the given grid", py::arg("raw_data"),
             py::arg("grid"), py::arg("max_peaks") = 0,
             py::arg("max_threads") = std::thread::hardware_concurrency())
        .def("warp_peaks", &PythonAPI::warp_peaks,
             "Warp peak lists to maximize the similarity with the given "
             "reference",
             py::arg("all_peaks"), py::arg("reference_index"), py::arg("slack"),
             py::arg("window_size"), py::arg("num_points"),
             py::arg("rt_expand_factor"), py::arg("peaks_per_window"))
        .def("find_similarity", &PythonAPI::find_similarity,
             "Find the similarity between two peak lists",
             py::arg("peak_list_a"), py::arg("peak_list_b"), py::arg("n_peaks"))
        .def("write_peaks", &PythonAPI::write_peaks,
             "Write the peaks to disk in a binary format", py::arg("peaks"),
             py::arg("file_name"))
        .def("write_metamatch_clusters", &PythonAPI::write_metamatch_clusters,
             "Write the metamatch_clusters to disk in a binary format",
             py::arg("metamatch_clusters"), py::arg("file_name"))
        .def("write_metamatch_peaks", &PythonAPI::write_metamatch_peaks,
             "Write the metamatch_peaks to disk in a binary format",
             py::arg("metamatch_peaks"), py::arg("file_name"))
        .def("write_linked_msms", &PythonAPI::write_linked_msms,
             "Write the linked_msms to disk in a binary format",
             py::arg("linked_msms"), py::arg("file_name"))
        .def("read_peaks", &PythonAPI::read_peaks,
             "Read the peaks from the binary peaks file", py::arg("file_name"))
        .def("read_metamatch_clusters", &PythonAPI::read_metamatch_clusters,
             "Read the metamatch_clusters from the binary metamatch_clusters "
             "file",
             py::arg("file_name"))
        .def("read_metamatch_peaks", &PythonAPI::read_metamatch_peaks,
             "Read the metamatch_peaks from the binary metamatch_peaks file",
             py::arg("file_name"))
        .def("read_raw_data", &PythonAPI::read_raw_data,
             "Read the raw_data from the binary raw_data file",
             py::arg("file_name"))
        .def("read_grid", &PythonAPI::read_grid,
             "Read the grid from the binary grid file", py::arg("file_name"))
        .def("read_linked_msms", &PythonAPI::read_linked_msms,
             "Read the linked_msms from the binary linked_msms file",
             py::arg("file_name"))
        .def("theoretical_isotopes_peptide",
             &FeatureDetection::theoretical_isotopes_peptide,
             "Calculate the theoretical isotopic distribution of a peptide",
             py::arg("sequence"), py::arg("charge_state"),
             py::arg("min_perc") = 0.01)
        .def("read_mzidentml", &PythonAPI::read_mzidentml,
             "Read identification data from the given mzIdentML file ",
             py::arg("file_name"))
        .def("read_ident_data", &PythonAPI::read_ident_data,
             "Read the ident_data from the binary ident_data file",
             py::arg("file_name"))
        .def("write_ident_data", &PythonAPI::write_ident_data,
             "Write the ident_data to disk in a binary format",
             py::arg("ident_data"), py::arg("file_name"))
        .def("read_features", &PythonAPI::read_features,
             "Read the feature from the binary feature file",
             py::arg("file_name"))
        .def("write_features", &PythonAPI::write_features,
             "Write the feature to disk in a binary format", py::arg("feature"),
             py::arg("file_name"))
        .def("read_feature_clusters", &PythonAPI::read_feature_clusters,
             "Read the feature_clusters from the binary feature_clusters file",
             py::arg("file_name"))
        .def("write_feature_clusters", &PythonAPI::write_feature_clusters,
             "Write the feature_clusters to disk in a binary format",
             py::arg("feature_clusters"), py::arg("file_name"))
        .def("perform_metamatch", &PythonAPI::perform_metamatch,
             "Perform metamatch for peak matching", py::arg("input"),
             py::arg("radius_mz"), py::arg("radius_rt"), py::arg("fraction"))
        .def("find_feature_clusters", &PythonAPI::find_feature_clusters,
             "Perform metamatch for feature matching", py::arg("group_ids"),
             py::arg("peaks"), py::arg("features"))
        .def("link_peaks", &Link::link_peaks, "Link msms events to peak ids",
             py::arg("peaks"), py::arg("raw_data"))
        .def("link_idents", &Link::link_idents,
             "Link msms events to spectrum identifications",
             py::arg("ident_data"), py::arg("raw_data"))
        .def("xic", &PythonAPI::xic, py::arg("raw_data"), py::arg("min_mz"),
             py::arg("max_mz"), py::arg("min_rt"), py::arg("max_rt"),
             py::arg("method") = "sum")
        .def("debug", &PythonAPI::debug, py::arg("ident_data"))
        .def("feature_detection", &FeatureDetection::feature_detection,
             "Link peaks as features", py::arg("peaks"),
             py::arg("raw_data_ms2"), py::arg("ident_data"),
             py::arg("link_table_idents"), py::arg("link_table_msms"),
             py::arg("discrepancy_threshold") = 0.25);
}
