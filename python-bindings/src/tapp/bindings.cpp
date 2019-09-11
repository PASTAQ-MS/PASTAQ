#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <thread>
#include <tuple>

#include "MIDAs/MIDAs.h"
#include "centroid/centroid.hpp"
#include "centroid/centroid_files.hpp"
#include "centroid/centroid_runners.hpp"
#include "grid/grid.hpp"
#include "grid/grid_files.hpp"
#include "grid/grid_serialize.hpp"
#include "grid/raw_data.hpp"
#include "grid/raw_data_serialize.hpp"
#include "grid/xml_reader.hpp"
#include "metamatch/metamatch.hpp"
#include "metamatch/metamatch_serialize.hpp"
#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "utils/search.hpp"
#include "utils/serialization.hpp"
#include "warp2d/warp2d.hpp"
#include "warp2d/warp2d_runners.hpp"

namespace py = pybind11;

namespace PythonAPI {
RawData::RawData read_mzxml(std::string &file_name, double min_mz,
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
    auto instrument_type = Instrument::Type::UNKNOWN;
    for (auto &ch : instrument_type_str) {
        ch = std::tolower(ch);
    }
    if (instrument_type_str == "orbitrap") {
        instrument_type = Instrument::Type::ORBITRAP;
    } else if (instrument_type_str == "tof") {
        instrument_type = Instrument::Type::TOF;
    } else if (instrument_type_str == "quad") {
        instrument_type = Instrument::Type::QUAD;
    } else if (instrument_type_str == "fticr") {
        instrument_type = Instrument::Type::FTICR;
    } else {
        std::ostringstream error_stream;
        error_stream << "the given instrument is not supported";
        throw std::invalid_argument(error_stream.str());
    }
    // Parse the polarity.
    auto polarity = RawData::Polarity::BOTH;
    for (auto &ch : polarity_str) {
        ch = std::tolower(ch);
    }
    if (polarity_str == "" || polarity_str == "both" || polarity_str == "+-" ||
        polarity_str == "-+") {
        polarity = RawData::Polarity::BOTH;
    } else if (polarity_str == "+" || polarity_str == "pos" ||
               polarity_str == "positive") {
        polarity = RawData::Polarity::POSITIVE;
    } else if (polarity_str == "-" || polarity_str == "neg" ||
               polarity_str == "negative") {
        polarity = RawData::Polarity::NEGATIVE;
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

    std::filesystem::path input_file = file_name;

    // Check for proper file extension.
    std::string extension = input_file.extension();
    std::string lowercase_extension = extension;
    for (auto &ch : lowercase_extension) {
        ch = std::tolower(ch);
    }
    if (lowercase_extension != ".mzxml") {
        std::ostringstream error_stream;
        error_stream << "invalid file type: expected 'mzXML' but given '"
                     << extension << "'";
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

std::string to_string(const Instrument::Type &instrument_type) {
    switch (instrument_type) {
        case Instrument::Type::QUAD:
            return "QUAD";
        case Instrument::Type::TOF:
            return "TOF";
        case Instrument::Type::FTICR:
            return "FTICR";
        case Instrument::Type::ORBITRAP:
            return "ORBITRAP";
        case Instrument::Type::UNKNOWN:
            return "UNKNOWN";
    };
    return "UNKNOWN";
}

std::string to_string(const RawData::Polarity &polarity) {
    switch (polarity) {
        case RawData::Polarity::POSITIVE:
            return "POSITIVE";
        case RawData::Polarity::NEGATIVE:
            return "NEGATIVE";
        case RawData::Polarity::BOTH:
            return "BOTH";
    };
    return "UNKNOWN";
}

struct LinkedMsms {
    size_t peak_id;
    size_t msms_id;
    size_t scan_index;
    double distance;
};

// NOTE: This algorithm relies on the peak vector to be sorted by id/height.
// TODO(alex): Move linkedmsms and serialization functions to separate
// namespace.
std::vector<LinkedMsms> link_msms(const std::vector<Centroid::Peak> &peaks,
                                  const RawData::RawData &raw_data) {
    // Index the peak list by m/z.
    struct PeakIndex {
        size_t id;
        double mz;
    };
    auto indices = std::vector<PeakIndex>(peaks.size());
    for (size_t i = 0; i < peaks.size(); ++i) {
        indices[i] = {peaks[i].id, peaks[i].local_max_mz};
    }

    // Sort mz index by mz.
    std::stable_sort(indices.begin(), indices.end(),
                     [](const PeakIndex &a, const PeakIndex &b) -> bool {
                         return a.mz < b.mz;
                     });

    // Flatten index into separate arrays.
    auto index_ids = std::vector<size_t>(indices.size());
    auto index_mzs = std::vector<double>(indices.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        index_ids[i] = indices[i].id;
        index_mzs[i] = indices[i].mz;
    }
    indices.clear();

    // Perform linkage of ms/ms events to closest peak.
    std::vector<LinkedMsms> link_table;
    for (size_t k = 0; k < raw_data.scans.size(); ++k) {
        const auto &scan = raw_data.scans[k];
        if (scan.ms_level != 2) {
            continue;
        }
        size_t event_id = scan.scan_number;
        double event_mz = scan.precursor_information.mz;
        double event_rt = scan.retention_time;
        double theoretical_sigma_mz = RawData::fwhm_to_sigma(
            RawData::theoretical_fwhm(raw_data, event_mz));

        // Find min_mz and loop until we reach the max_mz.
        double min_mz = event_mz - 3 * theoretical_sigma_mz;
        double max_mz = event_mz + 3 * theoretical_sigma_mz;
        size_t min_j = Search::lower_bound(index_mzs, min_mz);
        double min_distance = std::numeric_limits<double>::infinity();
        size_t peak_id = 0;
        for (size_t j = min_j; j < index_ids.size(); ++j) {
            size_t i = index_ids[j];
            if (peaks[i].local_max_mz > max_mz) {
                break;
            }
            double a = event_mz - peaks[i].local_max_mz;
            double b = event_rt - peaks[i].local_max_rt;
            double distance = std::sqrt(a * a + b * b);
            if (distance < min_distance) {
                min_distance = distance;
                peak_id = i;
            }
        }

        // Check if linked event is within 10 sigma of the minimum distance
        // peak.
        double roi_min_mz =
            peaks[peak_id].local_max_mz - 10 * peaks[peak_id].raw_roi_sigma_mz;
        double roi_max_mz =
            peaks[peak_id].local_max_mz + 10 * peaks[peak_id].raw_roi_sigma_mz;
        double roi_min_rt =
            peaks[peak_id].local_max_rt - 10 * peaks[peak_id].raw_roi_sigma_rt;
        double roi_max_rt =
            peaks[peak_id].local_max_rt + 10 * peaks[peak_id].raw_roi_sigma_rt;
        if (event_mz < roi_min_mz || event_mz > roi_max_mz ||
            event_rt < roi_min_rt || event_rt > roi_max_rt) {
            continue;
        }

        if (min_distance != std::numeric_limits<double>::infinity()) {
            link_table.push_back({peak_id, event_id, k, min_distance});
        }
    }

    // Sort link_table by peak_id.
    std::stable_sort(
        link_table.begin(), link_table.end(),
        [](const LinkedMsms &a, const LinkedMsms &b) -> bool {
            return (a.peak_id < b.peak_id) ||
                   ((a.peak_id == b.peak_id) && (a.distance < b.distance));
        });
    return link_table;
}

bool _read_linked_msms(std::istream &stream, LinkedMsms *link) {
    Serialization::read_uint64(stream, &link->peak_id);
    Serialization::read_uint64(stream, &link->msms_id);
    Serialization::read_uint64(stream, &link->scan_index);
    Serialization::read_double(stream, &link->distance);
    return stream.good();
}

bool _write_linked_msms(std::ostream &stream, const LinkedMsms &link) {
    Serialization::write_uint64(stream, link.peak_id);
    Serialization::write_uint64(stream, link.msms_id);
    Serialization::write_uint64(stream, link.scan_index);
    Serialization::write_double(stream, link.distance);
    return stream.good();
}

bool _read_linked_msms_table(std::istream &stream,
                             std::vector<LinkedMsms> *links) {
    uint64_t num_rows = 0;
    Serialization::read_uint64(stream, &num_rows);
    *links = std::vector<LinkedMsms>(num_rows);
    for (size_t i = 0; i < num_rows; ++i) {
        _read_linked_msms(stream, &(*links)[i]);
    }
    return stream.good();
}

bool _write_linked_msms_table(std::ostream &stream,
                              const std::vector<LinkedMsms> &links) {
    uint64_t num_rows = links.size();
    Serialization::write_uint64(stream, num_rows);
    for (size_t i = 0; i < num_rows; ++i) {
        _write_linked_msms(stream, links[i]);
    }
    return stream.good();
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
            Warp2D::Runners::Parallel::run(reference_peaks, peaks, parameters,
                                           std::thread::hardware_concurrency());
        all_warped_peaks.push_back(warped_peaks);
    }
    return all_warped_peaks;
}

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
    std::stable_sort(peak_list_a.begin(), peak_list_a.end(), sort_peaks);
    std::stable_sort(peak_list_b.begin(), peak_list_b.end(), sort_peaks);
    if (peak_list_a.size() > n_peaks) {
        peak_list_a.resize(n_peaks);
    }
    if (peak_list_b.size() > n_peaks) {
        peak_list_b.resize(n_peaks);
    }
    SimilarityResults results = {};
    results.self_a = Warp2D::similarity_2D(peak_list_a, peak_list_a);
    results.self_b = Warp2D::similarity_2D(peak_list_b, peak_list_b);
    results.overlap = Warp2D::similarity_2D(peak_list_a, peak_list_b);
    // Overlap / (GeometricMean(self_a, self_b))
    results.geometric_ratio =
        results.overlap / std::sqrt(results.self_a * results.self_b);
    // Harmonic mean of the ratios between
    // self_similarity/overlap_similarity
    results.mean_ratio =
        2 * results.overlap / (results.self_a + results.self_b);
    return results;
}

void to_csv(const std::vector<Centroid::Peak> &peaks, std::string &file_name) {
    std::filesystem::path output_file = file_name;

    // Open file stream.
    std::ofstream stream;
    stream.open(output_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open output file" << output_file;
        throw std::invalid_argument(error_stream.str());
    }

    if (!Centroid::Files::Csv::write_peaks(stream, peaks)) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't write the peaks into the output file"
                     << output_file;
        throw std::invalid_argument(error_stream.str());
    }
}

void write_raw_data(const RawData::RawData &raw_data, std::string &file_name) {
    std::filesystem::path output_file = file_name;

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

RawData::RawData read_raw_data(std::string &file_name) {
    std::filesystem::path input_file = file_name;

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

void write_mesh(const Grid::Mesh &mesh, std::string &file_name) {
    std::filesystem::path output_file = file_name;

    // Open file stream.
    std::ofstream stream;
    stream.open(output_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open output file" << output_file;
        throw std::invalid_argument(error_stream.str());
    }

    if (!Grid::Serialize::write_mesh(stream, mesh)) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't write the mesh into the output file"
                     << output_file;
        throw std::invalid_argument(error_stream.str());
    }
}

Grid::Mesh read_mesh(std::string &file_name) {
    std::filesystem::path input_file = file_name;

    // Open file stream.
    std::ifstream stream;
    stream.open(input_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }

    Grid::Mesh mesh;
    if (!Grid::Serialize::read_mesh(stream, &mesh)) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't write the mesh into the input file"
                     << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    return mesh;
}

void write_peaks(const std::vector<Centroid::Peak> &peaks,
                 std::string &file_name) {
    std::filesystem::path output_file = file_name;

    // Open file stream.
    std::ofstream stream;
    stream.open(output_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open output file" << output_file;
        throw std::invalid_argument(error_stream.str());
    }

    if (!Centroid::Files::Bpks::write_peaks(stream, peaks)) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't write the peaks into the output file"
                     << output_file;
        throw std::invalid_argument(error_stream.str());
    }
}

std::vector<Centroid::Peak> read_peaks(std::string &file_name) {
    std::filesystem::path input_file = file_name;

    // Open file stream.
    std::ifstream stream;
    stream.open(input_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }

    std::vector<Centroid::Peak> peaks;
    if (!Centroid::Files::Bpks::read_peaks(stream, &peaks)) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't write the peaks into the input file"
                     << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    return peaks;
}

void write_ident_data(const IdentData::IdentData &ident_data,
                      std::string &file_name) {
    std::filesystem::path output_file = file_name;

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

IdentData::IdentData read_ident_data(std::string &file_name) {
    std::filesystem::path input_file = file_name;

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
    std::string &file_name) {
    std::filesystem::path output_file = file_name;

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
    std::string &file_name) {
    std::filesystem::path input_file = file_name;

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
                           std::string &file_name) {
    std::filesystem::path output_file = file_name;

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

std::vector<MetaMatch::Peak> read_metamatch_peaks(std::string &file_name) {
    std::filesystem::path input_file = file_name;

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

void write_linked_msms(const std::vector<LinkedMsms> &linked_msms,
                       std::string &file_name) {
    std::filesystem::path output_file = file_name;

    // Open file stream.
    std::ofstream stream;
    stream.open(output_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open output file" << output_file;
        throw std::invalid_argument(error_stream.str());
    }

    if (!_write_linked_msms_table(stream, linked_msms)) {
        std::ostringstream error_stream;
        error_stream
            << "error: couldn't write the linked_msms into the output file"
            << output_file;
        throw std::invalid_argument(error_stream.str());
    }
}

std::vector<LinkedMsms> read_linked_msms(std::string &file_name) {
    std::filesystem::path input_file = file_name;

    // Open file stream.
    std::ifstream stream;
    stream.open(input_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }

    std::vector<LinkedMsms> linked_msms;
    if (!_read_linked_msms_table(stream, &linked_msms)) {
        std::ostringstream error_stream;
        error_stream
            << "error: couldn't write the linked_msms into the input file"
            << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    return linked_msms;
}

std::tuple<std::vector<double>, std::vector<double>>
theoretical_isotopes_peptide(std::string sequence, int8_t charge_state,
                             double min_perc) {
    auto midas = MIDAs(charge_state = charge_state);
    // NOTE: C00: Unmodified cysteine.
    //       C31: Carboxymethylation or Iodoacetic acid.
    //       C32: Carbamidomethylation or Iodoacetamide.
    //       C33: Pyridylethylation.
    midas.Initialize_Elemental_Composition(sequence, "C00", "H", "OH", 1);
    auto isotopes = midas.Coarse_Grained_Isotopic_Distribution();
    auto full_mzs = std::vector<double>(isotopes.size());
    auto probs = std::vector<double>(isotopes.size());
    double max_prob = 0;
    for (size_t i = 0; i < isotopes.size(); ++i) {
        full_mzs[i] = isotopes[i].mw / charge_state;
        probs[i] = isotopes[i].prob;
        if (probs[i] > max_prob) {
            max_prob = probs[i];
        }
    }
    // Normalize the probabilities to 0-1.
    std::vector<double> mzs;
    std::vector<double> perc;
    for (size_t i = 0; i < isotopes.size(); ++i) {
        probs[i] = probs[i] / max_prob;
        if (probs[i] > min_perc) {
            mzs.push_back(full_mzs[i]);
            perc.push_back(probs[i]);
        }
    }
    return {mzs, perc};
}

IdentData::IdentData read_mzidentml(std::string &file_name) {
    std::filesystem::path input_file = file_name;

    // Check for proper file extension.
    std::string extension = input_file.extension();
    std::string lowercase_extension = extension;
    for (auto &ch : lowercase_extension) {
        ch = std::tolower(ch);
    }
    if (lowercase_extension != ".mzid") {
        std::ostringstream error_stream;
        error_stream << "invalid file type: expected 'mzid' but given '"
                     << extension << "'";
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
    return XmlReader::read_mzidentml(stream);
}

struct Isotope {
    size_t id;
    double mz;
    double rt;
    double height;
    double intensity;
    double normalized_height;
    double expected_normalized_height;
};

struct LinkedPeptide {
    std::string sequence;
    std::string psm_id;
    size_t charge_state;
    double ident_rt;
    double ident_mz;
    std::vector<double> theoretical_isotopes_mz;
    std::vector<double> theoretical_isotopes_perc;

    std::vector<Isotope> linked_isotopes;
    double monoisotopic_height;
    double monoisotopic_intensity;
    double total_height;
    double total_intensity;
    double weighted_error;
};

std::vector<LinkedPeptide> link_identified_peptides(
    const std::vector<Centroid::Peak> &peaks,
    const std::vector<IdentData::SpectrumId> &identifications,
    double tolerance_rt, double minimum_isotope_perc) {
    std::vector<LinkedPeptide> linked_peptides;
    // Make a copy of the peaks and sort them by retention time.
    auto sorted_peaks = std::vector<Centroid::Peak>(peaks.size());
    for (size_t i = 0; i < peaks.size(); ++i) {
        sorted_peaks[i] = peaks[i];
    }
    std::stable_sort(sorted_peaks.begin(), sorted_peaks.end(),
                     [](const Centroid::Peak &p1, const Centroid::Peak &p2) {
                         return (p1.local_max_rt < p2.local_max_rt);
                     });

    // Extract sorted arrays for quick searching.
    auto sorted_rts = std::vector<double>(sorted_peaks.size());
    for (size_t i = 0; i < sorted_peaks.size(); ++i) {
        sorted_rts[i] = sorted_peaks[i].local_max_rt;
    }

    for (const auto &ident : identifications) {
        // Generate the theoretical isotope distribution for identified
        // peptides.
        auto [mzs, perc] = theoretical_isotopes_peptide(
            ident.sequence, ident.charge_state, minimum_isotope_perc);

        // Find the peaks within the rt tolerance range.
        double min_rt = ident.retention_time - tolerance_rt;
        double max_rt = ident.retention_time + tolerance_rt;
        size_t min_j = Search::lower_bound(sorted_rts, min_rt);
        size_t max_j = sorted_rts.size();
        std::vector<Centroid::Peak> peaks_in_range;
        for (size_t j = min_j; j < max_j; ++j) {
            if (sorted_rts[j] > max_rt) {
                break;
            }
            if (((sorted_peaks[j].local_max_mz +
                  sorted_peaks[j].raw_roi_sigma_mz) < mzs[0]) ||
                ((sorted_peaks[j].local_max_mz -
                  sorted_peaks[j].raw_roi_sigma_mz) > mzs[mzs.size() - 1])) {
                continue;
            }
            peaks_in_range.push_back(sorted_peaks[j]);
        }
        if (peaks_in_range.empty()) {
            continue;
        }

        // Sort the peaks in range by mz for easier searching.
        std::stable_sort(
            peaks_in_range.begin(), peaks_in_range.end(),
            [](const Centroid::Peak &p1, const Centroid::Peak &p2) {
                return (p1.local_max_mz < p2.local_max_mz);
            });

        // The reference node is the theoretial max relative to the
        // distribution.
        size_t reference_node_index = 0;

        // Create a graph the the potential isotope associations.
        std::vector<std::vector<Isotope>> isotopes_graph;
        for (size_t k = 0; k < mzs.size(); ++k) {
            std::vector<Isotope> candidates;
            double theoretical_mz = mzs[k];
            double theoretical_percentage = perc[k];
            if (theoretical_percentage == 1.0) {
                reference_node_index = k;
            }
            for (const auto &peak : peaks_in_range) {
                if ((peak.local_max_mz + peak.raw_roi_sigma_mz) <
                    theoretical_mz) {
                    continue;
                }
                if ((peak.local_max_mz - peak.raw_roi_sigma_mz) >
                    theoretical_mz) {
                    break;
                }
                candidates.push_back({
                    peak.id,
                    peak.local_max_mz,
                    peak.local_max_rt,
                    peak.local_max_height,
                    peak.raw_roi_total_intensity,
                    0.0,
                    0.0,
                });
            }
            isotopes_graph.push_back(candidates);
        }
        if (isotopes_graph.empty()) {
            // FIXME: Should we return NaNs instead?
            continue;
        }

        // In case more than one peak is linked to the reference mz isotope,
        // the sequence with the less matching error should be selected. In
        // order to do so, the heights for each candidate must be normalized
        // by the reference isotope height.
        std::vector<Isotope> selected_candidates;
        double weighted_error = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < isotopes_graph[reference_node_index].size();
             ++i) {
            const auto &ref_candidate = isotopes_graph[reference_node_index][i];
            std::vector<Isotope> normalized_candidates;
            for (size_t k = 0; k < mzs.size(); ++k) {
                if (k == reference_node_index) {
                    normalized_candidates.push_back(
                        {ref_candidate.id, ref_candidate.mz, ref_candidate.rt,
                         ref_candidate.height, ref_candidate.intensity, 1.0,
                         1.0});
                    continue;
                }

                // Find the best matching candidate for the selected
                // reference.
                double theoretical_percentage = perc[k];
                auto best_matching_candidate =
                    Isotope{0,
                            0.0,
                            0.0,
                            0.0,
                            0.0,
                            -std::numeric_limits<double>::infinity(),
                            theoretical_percentage};
                for (auto &candidate : isotopes_graph[k]) {
                    double normalized_height =
                        candidate.height / ref_candidate.height;
                    if (std::abs(candidate.normalized_height -
                                 theoretical_percentage) <
                        (std::abs(best_matching_candidate.normalized_height -
                                  theoretical_percentage))) {
                        best_matching_candidate =
                            Isotope{candidate.id,          candidate.mz,
                                    candidate.rt,          candidate.height,
                                    candidate.intensity,   normalized_height,
                                    theoretical_percentage};
                    }
                }
                normalized_candidates.push_back(best_matching_candidate);
            }
            // The weighted error for each normalized_candidate path can now
            // be evaluated. NOTE: A potential alternative to an error
            // function would be the Euclidean distance between the expected
            // theoretical peak and candidates for that node.
            double err_sum = 0.0;
            for (const auto &candidate : normalized_candidates) {
                err_sum += candidate.expected_normalized_height *
                           std::pow((candidate.expected_normalized_height -
                                     candidate.normalized_height),
                                    2);
            }
            if (err_sum < weighted_error) {
                weighted_error = err_sum;
                selected_candidates = normalized_candidates;
            }
        }
        if (selected_candidates.empty()) {
            continue;
        }
        LinkedPeptide linked_peptide = {};
        linked_peptide.sequence = ident.sequence;
        linked_peptide.psm_id = ident.id;
        linked_peptide.charge_state = ident.charge_state;
        linked_peptide.ident_rt = ident.retention_time;
        linked_peptide.ident_mz = ident.experimental_mz;
        linked_peptide.theoretical_isotopes_mz = mzs;
        linked_peptide.theoretical_isotopes_perc = perc;
        linked_peptide.linked_isotopes = selected_candidates;
        linked_peptide.monoisotopic_height = selected_candidates[0].height;
        linked_peptide.monoisotopic_intensity =
            selected_candidates[0].intensity;
        linked_peptide.weighted_error = weighted_error;
        for (const auto &candidate : selected_candidates) {
            linked_peptide.total_height += candidate.height;
            linked_peptide.total_intensity += candidate.intensity;
        }
        linked_peptides.push_back(linked_peptide);
    }
    return linked_peptides;
}

struct MetaMatchResults {
    std::vector<MetaMatch::Cluster> clusters;
    std::vector<MetaMatch::Peak> orphans;
};

MetaMatchResults perform_metamatch(
    // NOTE: [(class_0, peaks_0),...(class_i, peaks_i)]
    std::vector<std::tuple<size_t, std::vector<Centroid::Peak>>> input,
    double radius_mz, double radius_rt, double fraction) {
    MetaMatchResults results;

    // Create the ClassMaps.
    std::vector<MetaMatch::ClassMap> class_maps;
    std::vector<MetaMatch::Peak> metapeaks;
    size_t file_id = 0;
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

    py::class_<RawData::Polarity>(m, "Polarity")
        .def("__repr__", [](const RawData::Polarity &polarity) {
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
        .def("__repr__",
             [](const RawData::RawData &rd) {
                 return "RawData:\n> instrument_type: " +
                        PythonAPI::to_string(rd.instrument_type) +
                        "\n> resolution_ms1: " +
                        std::to_string(rd.resolution_ms1) +
                        "\n> resolution_msn: " +
                        std::to_string(rd.resolution_msn) +
                        "\n> reference_mz: " + std::to_string(rd.reference_mz) +
                        "\n> min_mz: " + std::to_string(rd.min_mz) +
                        "\n> max_mz: " + std::to_string(rd.max_mz) +
                        "\n> min_rt: " + std::to_string(rd.min_rt) +
                        "\n> max_rt: " + std::to_string(rd.max_rt) +
                        "\n> number of scans: " +
                        std::to_string(rd.scans.size());
             })
        .def("xic", &RawData::RawData::xic, py::arg("min_mz"),
             py::arg("max_mz"), py::arg("min_rt"), py::arg("max_rt"),
             py::arg("method") = "sum");

    py::class_<Grid::Mesh>(m, "Mesh")
        .def_readonly("n", &Grid::Mesh::n)
        .def_readonly("m", &Grid::Mesh::m)
        .def_readonly("data", &Grid::Mesh::data)
        .def_readonly("bins_mz", &Grid::Mesh::bins_mz)
        .def_readonly("bins_rt", &Grid::Mesh::bins_rt)
        .def("dump", &PythonAPI::write_mesh)
        .def("__repr__", [](const Grid::Mesh &s) {
            return "Mesh <n: " + std::to_string(s.n) +
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

    py::class_<Centroid::Peak>(m, "Peak")
        .def_readonly("id", &Centroid::Peak::id)
        .def_readonly("local_max_mz", &Centroid::Peak::local_max_mz)
        .def_readonly("local_max_rt", &Centroid::Peak::local_max_rt)
        .def_readonly("local_max_height", &Centroid::Peak::local_max_height)
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
        .def("xic", &Centroid::Peak::xic, py::arg("raw_data"),
             py::arg("method") = "sum")
        .def("__repr__", [](const Centroid::Peak &p) {
            return "Peak <id: " + std::to_string(p.id) +
                   ", local_max_mz: " + std::to_string(p.local_max_mz) +
                   ", local_max_rt: " + std::to_string(p.local_max_rt) +
                   ", local_max_height: " + std::to_string(p.local_max_height) +
                   ", raw_roi_sigma_mz: " + std::to_string(p.raw_roi_sigma_mz) +
                   ", raw_roi_sigma_rt: " + std::to_string(p.raw_roi_sigma_rt) +
                   ", raw_roi_num_points: " +
                   std::to_string(p.raw_roi_num_points) +
                   ", raw_roi_num_scans: " +
                   std::to_string(p.raw_roi_num_scans) + ">";
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

    py::class_<PythonAPI::LinkedPeptide>(m, "LinkedPeptide")
        .def_readonly("sequence", &PythonAPI::LinkedPeptide::sequence)
        .def_readonly("psm_id", &PythonAPI::LinkedPeptide::psm_id)
        .def_readonly("charge_state", &PythonAPI::LinkedPeptide::charge_state)
        .def_readonly("ident_rt", &PythonAPI::LinkedPeptide::ident_rt)
        .def_readonly("ident_mz", &PythonAPI::LinkedPeptide::ident_mz)
        .def_readonly("theoretical_isotopes_mz",
                      &PythonAPI::LinkedPeptide::theoretical_isotopes_mz)
        .def_readonly("theoretical_isotopes_perc",
                      &PythonAPI::LinkedPeptide::theoretical_isotopes_perc)
        .def_readonly("linked_isotopes",
                      &PythonAPI::LinkedPeptide::linked_isotopes)
        .def_readonly("monoisotopic_height",
                      &PythonAPI::LinkedPeptide::monoisotopic_height)
        .def_readonly("monoisotopic_intensity",
                      &PythonAPI::LinkedPeptide::monoisotopic_intensity)
        .def_readonly("total_height", &PythonAPI::LinkedPeptide::total_height)
        .def_readonly("total_intensity",
                      &PythonAPI::LinkedPeptide::total_intensity)
        .def_readonly("weighted_error",
                      &PythonAPI::LinkedPeptide::weighted_error)
        .def("__repr__", [](const PythonAPI::LinkedPeptide &s) {
            return s.sequence + "_" + std::to_string(s.charge_state);
        });

    py::class_<PythonAPI::Isotope>(m, "Isotope")
        .def_readonly("id", &PythonAPI::Isotope::id)
        .def_readonly("mz", &PythonAPI::Isotope::mz)
        .def_readonly("rt", &PythonAPI::Isotope::rt)
        .def_readonly("height", &PythonAPI::Isotope::height)
        .def_readonly("intensity", &PythonAPI::Isotope::intensity)
        .def_readonly("normalized_height",
                      &PythonAPI::Isotope::normalized_height)
        .def_readonly("expected_normalized_height",
                      &PythonAPI::Isotope::expected_normalized_height)
        .def("__repr__",
             [](const PythonAPI::Isotope &s) { return std::to_string(s.id); });

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

    py::class_<PythonAPI::LinkedMsms>(m, "LinkedMsms")
        .def_readonly("peak_id", &PythonAPI::LinkedMsms::peak_id)
        .def_readonly("msms_id", &PythonAPI::LinkedMsms::msms_id)
        .def_readonly("distance", &PythonAPI::LinkedMsms::distance)
        .def("__repr__", [](const PythonAPI::LinkedMsms &p) {
            return "LinkedMsms <peak_id: " + std::to_string(p.peak_id) +
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
        .def("resample", &Grid::resample,
             "Resample the raw data into a smoothed warped grid",
             py::arg("raw_data"), py::arg("num_mz") = 10,
             py::arg("num_rt") = 10, py::arg("smoothing_coef_mz") = 0.5,
             py::arg("smoothing_coef_rt") = 0.5)
        .def("find_raw_points", &RawData::find_raw_points,
             "Save the fitted peaks as a bpks file", py::arg("raw_data"),
             py::arg("min_mz"), py::arg("max_mz"), py::arg("min_rt"),
             py::arg("max_rt"))
        .def("find_peaks", &Centroid::Runners::Parallel::run,
             "Find all peaks in the given mesh", py::arg("raw_data"),
             py::arg("mesh"), py::arg("max_peaks") = 0)
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
        .def("to_csv", &PythonAPI::to_csv,
             "Write the peaks to disk in csv format (compatibility)",
             py::arg("peaks"), py::arg("file_name"))
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
        .def("read_mesh", &PythonAPI::read_mesh,
             "Read the mesh from the binary mesh file", py::arg("file_name"))
        .def("read_linked_msms", &PythonAPI::read_linked_msms,
             "Read the linked_msms from the binary linked_msms file",
             py::arg("file_name"))
        .def("theoretical_isotopes_peptide",
             &PythonAPI::theoretical_isotopes_peptide,
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
        .def("perform_metamatch", &PythonAPI::perform_metamatch,
             "Perform metamatch for peak matching", py::arg("input"),
             py::arg("radius_mz"), py::arg("radius_rt"), py::arg("fraction"))
        .def("link_msms", &PythonAPI::link_msms, "Link msms events to peak ids",
             py::arg("peaks"), py::arg("raw_data"))
        .def("link_identified_peptides", &PythonAPI::link_identified_peptides,
             "DEBUG", py::arg("peaks"), py::arg("identifications"),
             py::arg("tolerance_rt"), py::arg("minimum_isotope_perc"));
}
