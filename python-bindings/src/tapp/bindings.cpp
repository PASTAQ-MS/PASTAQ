#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <tuple>

#include "centroid/centroid.hpp"
#include "centroid/centroid_files.hpp"
#include "grid/grid.hpp"
#include "grid/grid_files.hpp"
#include "grid/xml_reader.hpp"
#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

namespace py = pybind11;

namespace PythonAPI {
RawData::RawData read_mzxml(std::string file_name, double min_mz, double max_mz,
                            double min_rt, double max_rt,
                            std::string instrument_type_str,
                            double resolution_ms1, double resolution_msn,
                            double reference_mz, double fwhm_rt,
                            std::string polarity_str) {
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
        resolution_msn, reference_mz, polarity);
    if (!raw_data) {
        std::ostringstream error_stream;
        error_stream << "error: an error occurred when reading the file"
                     << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    raw_data.value().fwhm_rt = fwhm_rt;

    return raw_data.value();
}

struct Roi {
    double min_mz;
    double max_mz;
    double min_rt;
    double max_rt;
};

struct RoiIndices {
    size_t min_j;
    size_t max_j;

    std::vector<size_t> min_i;
    std::vector<size_t> max_i;
};

struct RawPoint {
    double rt;
    double mz;
    double intensity;
};

struct RawPoints {
    std::vector<double> rt;
    std::vector<double> mz;
    std::vector<double> intensity;
    size_t num_points;
};

RawPoints find_raw_points(const RawData::RawData &raw_data, double min_mz,
                          double max_mz, double min_rt, double max_rt) {
    RawPoints raw_points;
    const auto &scans = raw_data.scans;
    if (scans.size() == 0) {
        std::ostringstream error_stream;
        error_stream << "the given raw_data is empty";
        throw std::invalid_argument(error_stream.str());
    }

    // Find scan indices.
    if (min_rt < raw_data.min_rt) {
        min_rt = raw_data.min_rt;
    }
    if (max_rt > raw_data.max_rt) {
        max_rt = raw_data.max_rt;
    }

    // Binary search for lower rt bound.
    size_t min_j = 0;
    size_t max_j = scans.size();
    size_t l = min_j;
    size_t r = max_j - 1;
    while (l <= r) {
        min_j = (l + r) / 2;
        if (scans[min_j].retention_time < min_rt) {
            l = min_j + 1;
        } else if (scans[min_j].retention_time > min_rt) {
            r = min_j - 1;
        } else {
            break;
        }
        if (min_j == 0) {
            break;
        }
    }
    for (size_t j = min_j; j < max_j; ++j) {
        const auto &scan = scans[j];
        if (scan.num_points == 0) {
            continue;
        }
        if (scan.retention_time > max_rt) {
            break;
        }
        // Binary search for lower mz bound.
        double internal_min_mz = min_mz;
        double internal_max_mz = max_mz;
        if (internal_min_mz < scan.mz[0]) {
            internal_min_mz = scan.mz[0];
        }
        if (internal_max_mz > scan.mz[scan.num_points - 1]) {
            internal_max_mz = scan.mz[scan.num_points - 1];
        }
        // Binary search for lower bound.
        size_t min_i = 0;
        size_t max_i = scan.num_points;
        size_t l = min_i;
        size_t r = max_i - 1;
        while (l <= r) {
            min_i = (l + r) / 2;
            if (scan.mz[min_i] < internal_min_mz) {
                l = min_i + 1;
            } else if (scan.mz[min_i] > internal_min_mz) {
                r = min_i - 1;
            } else {
                break;
            }
            if (min_i == 0) {
                break;
            }
        }
        for (size_t i = min_i; i < max_i; ++i) {
            if (scan.mz[i] > internal_max_mz) {
                break;
            }

            raw_points.rt.push_back(scan.retention_time);
            raw_points.mz.push_back(scan.mz[i]);
            raw_points.intensity.push_back(scan.intensity[i]);
            ++raw_points.num_points;
        }
    }

    if (raw_points.num_points == 0) {
        std::ostringstream error_stream;
        error_stream << "couldn't find raw_data points on the given ROI";
        throw std::invalid_argument(error_stream.str());
    }

    return raw_points;
}

uint64_t x_index(const RawData::RawData &raw_data, double mz, uint64_t k) {
    double fwhm_ref = raw_data.reference_mz / raw_data.resolution_ms1;
    switch (raw_data.instrument_type) {
        case Instrument::ORBITRAP: {
            double a = fwhm_ref / std::pow(raw_data.reference_mz, 1.5);
            double b = (1 / std::sqrt(raw_data.min_mz) - 1 / std::sqrt(mz));
            return static_cast<uint64_t>(k * 2 / a * b);
        } break;
        case Instrument::FTICR: {
            double a = 1 - raw_data.min_mz / mz;
            double b = raw_data.reference_mz * raw_data.reference_mz;
            double c = fwhm_ref * raw_data.min_mz;
            return static_cast<uint64_t>(k * a * b / c);
        } break;
        case Instrument::TOF: {
            return static_cast<uint64_t>(k * raw_data.reference_mz / fwhm_ref *
                                         std::log(mz / raw_data.min_mz));
        } break;
        case Instrument::QUAD: {
            // Same as the regular grid.
            return static_cast<uint64_t>(k * (mz - raw_data.min_mz) / fwhm_ref);
        } break;
        case Instrument::UNKNOWN: {
            assert(false);  // Can't handle unknown instruments.
        } break;
    }
    assert(false);  // Can't handle unknown instruments.
    return 0;
}

uint64_t y_index(const RawData::RawData &raw_data, double rt, uint64_t k) {
    double delta_rt = raw_data.fwhm_rt / k;
    return std::ceil((rt - raw_data.min_rt) / delta_rt);
}

std::tuple<uint64_t, uint64_t> calculate_dimensions(
    const RawData::RawData &raw_data, uint64_t k, uint64_t t) {
    return std::tuple<uint64_t, uint64_t>(
        x_index(raw_data, raw_data.max_mz, k) + 1,
        y_index(raw_data, raw_data.max_rt, t) + 1);
}

double fwhm_at(const RawData::RawData &raw_data, double mz) {
    double e = 0;
    switch (raw_data.instrument_type) {
        case Instrument::ORBITRAP: {
            e = 1.5;
        } break;
        case Instrument::FTICR: {
            e = 2;
        } break;
        case Instrument::TOF: {
            e = 1;
        } break;
        case Instrument::QUAD: {
            e = 0;
        } break;
        case Instrument::UNKNOWN: {
            assert(false);  // Can't handle unknown instruments.
        } break;
    }
    double mz_ref = raw_data.reference_mz;
    double fwhm_ref = mz_ref / raw_data.resolution_ms1;
    return fwhm_ref * std::pow(mz / mz_ref, e);
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

struct Mesh {
    uint64_t n;
    uint64_t m;
    std::vector<double> matrix;
    std::vector<double> bins_mz;
    std::vector<double> bins_rt;

    // Dumps the mesh into a binary .dat file.
    void save(std::string file_name) {
        // Open file stream.
        std::filesystem::path output_file = file_name;
        std::ofstream stream;
        stream.open(output_file);
        if (!stream) {
            std::ostringstream error_stream;
            error_stream << "error: couldn't open output file" << output_file;
            throw std::invalid_argument(error_stream.str());
        }

        // TODO: Error checking, out of bounds, correctness, etc.
        auto parameters = Grid::Parameters{};
        parameters.dimensions.n = n;
        parameters.dimensions.m = m;
        parameters.bounds.min_rt = bins_rt[0];
        parameters.bounds.max_rt = bins_rt[m - 1];
        parameters.bounds.min_mz = bins_mz[0];
        parameters.bounds.max_mz = bins_mz[n - 1];
        // FIXME: For now... Should we store instrument_type in Mesh?
        parameters.instrument_type = Instrument::Type::ORBITRAP;
        if (!Grid::Files::Dat::write(stream, matrix, parameters)) {
            std::cout << "error: the grid could not be saved properly"
                      << std::endl;
            return;
        }
    }
};

double fwhm_to_sigma(double fwhm) {
    return fwhm / (2 * std::sqrt(2 * std::log(2)));
}

// Applies a 2D kernel smoothing. The smoothing is performed in two passes.
// First the raw data points are mapped into a 2D matrix by splatting them into
// a matrix. Sparse areas might result in artifacts when the data is noisy, for
// this reason, the data is smoothed again.
//
// Since multiple passes of a Gaussian smoothing is equivalent to a single
// pass with `sigma = sqrt(2) * sigma_pass`, we adjust the sigmas for each pass
// accordingly.
Mesh resample(const RawData::RawData &raw_data, uint64_t num_samples_mz,
              uint64_t num_samples_rt, double smoothing_coef_mz,
              double smoothing_coef_rt) {
    // Initialize the Mesh object.
    auto [n, m] =
        calculate_dimensions(raw_data, num_samples_mz, num_samples_rt);
    Mesh mesh;
    mesh.n = n;
    mesh.m = m;
    mesh.matrix = std::vector<double>(n * m);
    auto weights = std::vector<double>(n * m);
    mesh.bins_mz = std::vector<double>(n);
    mesh.bins_rt = std::vector<double>(m);

    // Generate bins_mz.
    double mz_ref = raw_data.reference_mz;
    double fwhm_ref = raw_data.reference_mz / raw_data.resolution_ms1;
    for (size_t i = 0; i < n; ++i) {
        switch (raw_data.instrument_type) {
            case Instrument::ORBITRAP: {
                double a = 1 / std::sqrt(raw_data.min_mz);
                double b =
                    fwhm_ref / std::pow(mz_ref, 1.5) * i / 2 / num_samples_mz;
                double c = a - b;
                mesh.bins_mz[i] = 1 / (c * c);
            } break;
            case Instrument::FTICR: {
                double a = fwhm_ref * raw_data.min_mz;
                double b = mz_ref * mz_ref;
                mesh.bins_mz[i] =
                    raw_data.min_mz / (1 - (a / b) * i / num_samples_mz);
            } break;
            case Instrument::TOF: {
                mesh.bins_mz[i] =
                    raw_data.min_mz *
                    std::exp(fwhm_ref / mz_ref * i / num_samples_mz);
            } break;
            case Instrument::QUAD: {
                double delta_mz = (raw_data.max_mz - raw_data.min_mz) /
                                  static_cast<double>(mesh.n - 1);
                mesh.bins_mz[i] =
                    raw_data.min_mz + delta_mz * i / num_samples_mz;
            } break;
            case Instrument::UNKNOWN: {
                assert(false);  // Can't handle unknown instruments.
            } break;
        }
    }

    // Generate bins_rt.
    double delta_rt = (raw_data.max_rt - raw_data.min_rt) / (m - 1);
    for (size_t j = 0; j < m; ++j) {
        mesh.bins_rt[j] = raw_data.min_rt + delta_rt * j;
    }

    // Pre-calculate the smoothing sigma values for all bins of the grid.
    double sigma_rt =
        fwhm_to_sigma(raw_data.fwhm_rt) * smoothing_coef_rt / std::sqrt(2);
    auto sigma_mz_vec = std::vector<double>(n);
    for (size_t i = 0; i < n; ++i) {
        sigma_mz_vec[i] = fwhm_to_sigma(fwhm_at(raw_data, mesh.bins_mz[i])) *
                          smoothing_coef_mz / std::sqrt(2);
    }

    // Pre-calculate the kernel half widths for rt and all mzs.
    //
    // Since sigma_rt is constant, the size of the kernel will be the same
    // for the entire rt range.
    uint64_t rt_kernel_hw = 3 * sigma_rt / delta_rt;
    auto mz_kernel_hw = std::vector<uint64_t>(n);
    for (size_t i = 0; i < n; ++i) {
        double sigma_mz = sigma_mz_vec[i];
        double delta_mz = 0;
        if (i == 0) {
            delta_mz = mesh.bins_mz[i + 1] - mesh.bins_mz[i];
        } else {
            delta_mz = mesh.bins_mz[i] - mesh.bins_mz[i - 1];
        }
        mz_kernel_hw[i] = 3 * sigma_mz / delta_mz;
    }

    // Gaussian splatting.
    for (size_t i = 0; i < raw_data.scans.size(); ++i) {
        const auto &scan = raw_data.scans[i];
        double current_rt = scan.retention_time;

        // Find the bin for the current retention time.
        size_t index_rt = y_index(raw_data, current_rt, num_samples_rt);

        // Find the min/max indexes for the rt kernel.
        size_t j_min = 0;
        if (index_rt >= rt_kernel_hw) {
            j_min = index_rt - rt_kernel_hw;
        }
        size_t j_max = mesh.m - 1;
        if ((index_rt + rt_kernel_hw) < mesh.m) {
            j_max = index_rt + rt_kernel_hw;
        }

        for (size_t k = 0; k < scan.num_points; ++k) {
            double current_intensity = scan.intensity[k];
            double current_mz = scan.mz[k];

            // Find the bin for the current mz.
            size_t index_mz = x_index(raw_data, current_mz, num_samples_mz);

            double sigma_mz = sigma_mz_vec[index_mz];

            // Find the min/max indexes for the mz kernel.
            size_t i_min = 0;
            if (index_mz >= mz_kernel_hw[index_mz]) {
                i_min = index_mz - mz_kernel_hw[index_mz];
            }
            size_t i_max = mesh.n - 1;
            if ((index_mz + mz_kernel_hw[index_mz]) < mesh.n) {
                i_max = index_mz + mz_kernel_hw[index_mz];
            }

            for (size_t j = j_min; j <= j_max; ++j) {
                for (size_t i = i_min; i <= i_max; ++i) {
                    double x = mesh.bins_mz[i];
                    double y = mesh.bins_rt[j];

                    // Calculate the Gaussian weight for this point.
                    double a = (x - current_mz) / sigma_mz;
                    double b = (y - current_rt) / sigma_rt;
                    double weight = std::exp(-0.5 * (a * a + b * b));

                    mesh.matrix[i + j * n] += weight * current_intensity;
                    weights[i + j * n] += weight;
                }
            }
        }
    }
    for (size_t i = 0; i < (n * m); ++i) {
        double weight = weights[i];
        if (weight == 0) {
            weight = 1;
        }
        mesh.matrix[i] = mesh.matrix[i] / weight;
    }

    // Gaussian smoothing.
    //
    // The Gaussian 2D filter is separable. We obtain the same result with
    // faster performance by applying two 1D kernel convolutions instead. This
    // is specially noticeable on the full image.
    {
        auto smoothed_matrix = std::vector<double>(mesh.n * mesh.m);

        // Retention time smoothing.
        for (size_t j = 0; j < mesh.m; ++j) {
            double current_rt = mesh.bins_rt[j];
            size_t min_k = 0;
            if (j >= rt_kernel_hw) {
                min_k = j - rt_kernel_hw;
            }
            size_t max_k = mesh.m - 1;
            if ((j + rt_kernel_hw) < mesh.m) {
                max_k = j + rt_kernel_hw;
            }
            for (size_t i = 0; i < mesh.n; ++i) {
                double sum_weights = 0;
                double sum_weighted_values = 0;
                for (size_t k = min_k; k <= max_k; ++k) {
                    double a = (current_rt - mesh.bins_rt[k]) / sigma_rt;
                    double weight = std::exp(-0.5 * (a * a));
                    sum_weights += weight;
                    sum_weighted_values += weight * mesh.matrix[i + k * mesh.n];
                }
                smoothed_matrix[i + j * mesh.n] =
                    sum_weighted_values / sum_weights;
            }
        }
        mesh.matrix = smoothed_matrix;
    }
    {
        auto smoothed_matrix = std::vector<double>(mesh.n * mesh.m);

        // mz smoothing.
        //
        // Since sigma_mz is not constant, we need to calculate the kernel half
        // width for each potential mz value.
        for (size_t i = 0; i < mesh.n; ++i) {
            double sigma_mz = sigma_mz_vec[i];
            double current_mz = mesh.bins_mz[i];

            size_t min_k = 0;
            if (i >= mz_kernel_hw[i]) {
                min_k = i - mz_kernel_hw[i];
            }
            size_t max_k = mesh.n - 1;
            if ((i + mz_kernel_hw[i]) < mesh.n) {
                max_k = i + mz_kernel_hw[i];
            }
            for (size_t j = 0; j < mesh.m; ++j) {
                double sum_weights = 0;
                double sum_weighted_values = 0;
                for (size_t k = min_k; k <= max_k; ++k) {
                    double a = (current_mz - mesh.bins_mz[k]) / sigma_mz;
                    double weight = std::exp(-0.5 * (a * a));
                    sum_weights += weight;
                    sum_weighted_values += weight * mesh.matrix[k + j * mesh.n];
                }
                smoothed_matrix[i + j * mesh.n] =
                    sum_weighted_values / sum_weights;
            }
        }
        mesh.matrix = smoothed_matrix;
    }
    return mesh;
}

// TODO: For now returns a tuple, we should maybe return a definite struct.
//
//     returns: (i,j,mz,rt,intensity)
//
std::vector<std::tuple<uint64_t, uint64_t, double, double, double>>
find_local_max(const Mesh &mesh) {
    std::vector<std::tuple<uint64_t, uint64_t, double, double, double>> points;
    // FIXME: This is performed in O(n^2), but using the divide and conquer
    // strategy we might achieve O(n * log(n)) or lower.
    // FIXME: Also, we should consider the corner case where neighbours are
    // exactly equal, both should be considered a local maxima and the average
    // of mz and rt should be reported.
    for (size_t j = 1; j < mesh.m - 1; ++j) {
        for (size_t i = 1; i < mesh.n - 1; ++i) {
            int index = i + j * mesh.n;

            // NOTE(alex): The definition of a local maxima in a 2D space might
            // have different interpretations. i.e. We can select the 8
            // neighbours and the local maxima will be marked if all points are
            // below the central value. Alternatively, only a number N of
            // neighbours can be used, for example only the 4 cardinal
            // directions from the value under study.
            //
            // ----------------------------------------------
            // |              | top_value    |              |
            // ----------------------------------------------
            // | left_value   | value        | right_value  |
            // ----------------------------------------------
            // |              | bottom_value |              |
            // ----------------------------------------------
            double value = mesh.matrix[index];
            double right_value = mesh.matrix[index + 1];
            double left_value = mesh.matrix[index - 1];
            double top_value = mesh.matrix[index - mesh.n];
            double bottom_value = mesh.matrix[index + mesh.n];

            if ((value != 0) && (value > left_value) && (value > right_value) &&
                (value > top_value) && (value > bottom_value)) {
                points.push_back(
                    {i, j, mesh.bins_mz[i], mesh.bins_rt[j], value});
            }
        }
    }

    return points;
}

// FIXME: Terrible!
// Tuple:
//
//     0      ,  1 ,  2 ,  3  ,  4  ,  5        ,  6        ,  7
//     height ,  i ,  j ,  mz ,  rt ,  sigma_mz ,  sigma_rt ,  total_intensity
//
void save_fitted_peaks(
    const std::vector<std::tuple<double, double, double, double, double, double,
                                 double, double>> &fitted_peaks,
    std::string file_name) {
    std::vector<Centroid::Peak> peaks;
    for (const auto &fitted_peak : fitted_peaks) {
        Centroid::Peak peak = {};
        peak.height = std::get<0>(fitted_peak);
        peak.i = std::get<1>(fitted_peak);
        peak.j = std::get<2>(fitted_peak);
        peak.mz = std::get<3>(fitted_peak);
        peak.rt = std::get<4>(fitted_peak);
        peak.sigma_mz = std::get<5>(fitted_peak);
        peak.sigma_rt = std::get<6>(fitted_peak);
        peak.total_intensity = std::get<7>(fitted_peak);
        peak.mz_centroid = std::get<3>(fitted_peak);
        peak.rt_centroid = std::get<4>(fitted_peak);
        peak.height_centroid = std::get<0>(fitted_peak);
        peak.total_intensity_centroid = std::get<7>(fitted_peak);
        peaks.push_back(peak);
    }

    std::filesystem::path output_file = file_name;
    std::ofstream stream;
    stream.open(output_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open output file" << output_file;
        throw std::invalid_argument(error_stream.str());
    }
    Grid::Parameters grid_params;
    if (!Centroid::Files::Bpks::write_peaks(stream, grid_params, peaks)) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't peaks into file " << file_name;
        throw std::invalid_argument(error_stream.str());
    }

    std::cout << "Saving peaks to disk in csv..." << std::endl;
    auto csv_outfile_name = output_file.filename().replace_extension(".csv");
    std::ofstream csv_outfile_stream;
    csv_outfile_stream.open(csv_outfile_name, std::ios::out | std::ios::binary);
    std::cout << "Sorting peaks by height (centroid)..." << std::endl;
    auto sort_peaks = [](const Centroid::Peak &p1,
                         const Centroid::Peak &p2) -> bool {
        return (p1.height_centroid > p2.height_centroid) ||
               ((p1.height_centroid == p2.height_centroid) &&
                (p1.total_intensity_centroid > p2.total_intensity_centroid));
    };
    std::stable_sort(peaks.begin(), peaks.end(), sort_peaks);
    if (!Centroid::Files::Csv::write_peaks(csv_outfile_stream, peaks)) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't write peaks into file "
                     << csv_outfile_name;
        throw std::invalid_argument(error_stream.str());
    }
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
        .def_readonly("polarity", &RawData::Scan::polarity);

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
        .def_readonly("instrument_type", &RawData::RawData::instrument_type)
        .def_readonly("resolution_ms1", &RawData::RawData::resolution_ms1)
        .def_readonly("resolution_msn", &RawData::RawData::resolution_msn)
        .def_readonly("reference_mz", &RawData::RawData::reference_mz)
        .def_readonly("min_mz", &RawData::RawData::min_mz)
        .def_readonly("max_mz", &RawData::RawData::max_mz)
        .def_readonly("min_rt", &RawData::RawData::min_rt)
        .def_readonly("max_rt", &RawData::RawData::max_rt)
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

    py::class_<PythonAPI::Mesh>(m, "Mesh")
        .def_readonly("n", &PythonAPI::Mesh::n)
        .def_readonly("m", &PythonAPI::Mesh::m)
        .def_readonly("matrix", &PythonAPI::Mesh::matrix)
        .def_readonly("bins_mz", &PythonAPI::Mesh::bins_mz)
        .def_readonly("bins_rt", &PythonAPI::Mesh::bins_rt)
        .def("save", &PythonAPI::Mesh::save, py::arg("file_name"));

    py::class_<PythonAPI::RawPoints>(m, "RawPoints")
        .def_readonly("rt", &PythonAPI::RawPoints::rt)
        .def_readonly("mz", &PythonAPI::RawPoints::mz)
        .def_readonly("intensity", &PythonAPI::RawPoints::intensity);

    // Functions.
    m.def("read_mzxml", &PythonAPI::read_mzxml,
          "Read raw data from the given mzXML file ", py::arg("file_name"),
          py::arg("min_mz") = -1.0, py::arg("max_mz") = -1.0,
          py::arg("min_rt") = -1.0, py::arg("max_rt") = -1.0,
          py::arg("instrument_type") = "", py::arg("resolution_ms1"),
          py::arg("resolution_msn"), py::arg("reference_mz"),
          py::arg("fwhm_rt"), py::arg("polarity") = "")
        .def("calculate_dimensions", &PythonAPI::calculate_dimensions,
             "Calculate the grid parameters for the given raw file",
             py::arg("raw_data"), py::arg("num_mz") = 10,
             py::arg("num_rt") = 10)
        .def("fwhm_at", &PythonAPI::fwhm_at,
             "Calculate the width of the peak at the given m/z for the given "
             "raw file",
             py::arg("raw_data"), py::arg("mz"))
        .def("resample", &PythonAPI::resample,
             "Resample the raw data into a smoothed warped grid",
             py::arg("raw_data"), py::arg("num_mz") = 10,
             py::arg("num_rt") = 10, py::arg("smoothing_coef_mz") = 0.5,
             py::arg("smoothing_coef_rt") = 0.5)
        .def("find_local_max", &PythonAPI::find_local_max,
             "Find all local maxima in the given mesh", py::arg("mesh"))
        .def("save_fitted_peaks", &PythonAPI::save_fitted_peaks,
             "Save the fitted peaks as a bpks file", py::arg("fitted_peaks"),
             py::arg("file_name"))
        .def("find_raw_points", &PythonAPI::find_raw_points,
             "Save the fitted peaks as a bpks file", py::arg("raw_data"),
             py::arg("min_mz"), py::arg("max_mz"), py::arg("min_rt"),
             py::arg("max_rt"));
}
