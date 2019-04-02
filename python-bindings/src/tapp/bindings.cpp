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
#include "warp2d/warp2d.hpp"
#include "warp2d/warp2d_runners.hpp"

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
    uint64_t n;  // Number of mz sampling points.
    uint64_t m;  // Number of rt sampling points.
    uint64_t k;  // Number of sampling points per FWHM in mz.
    uint64_t t;  // Number of sampling points per FWHM in rt.
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
    mesh.k = num_samples_mz;
    mesh.t = num_samples_rt;
    mesh.matrix = std::vector<double>(n * m);
    mesh.bins_mz = std::vector<double>(n);
    mesh.bins_rt = std::vector<double>(m);

    // Generate bins_mz.
    double mz_ref = raw_data.reference_mz;
    double fwhm_ref = raw_data.reference_mz / raw_data.resolution_ms1;
    for (size_t i = 0; i < n; ++i) {
        switch (raw_data.instrument_type) {
            case Instrument::ORBITRAP: {
                double a = 1 / std::sqrt(raw_data.min_mz);
                double b = fwhm_ref / std::pow(mz_ref, 1.5) * i / 2 / mesh.k;
                double c = a - b;
                mesh.bins_mz[i] = 1 / (c * c);
            } break;
            case Instrument::FTICR: {
                double a = fwhm_ref * raw_data.min_mz;
                double b = mz_ref * mz_ref;
                mesh.bins_mz[i] = raw_data.min_mz / (1 - (a / b) * i / mesh.k);
            } break;
            case Instrument::TOF: {
                mesh.bins_mz[i] =
                    raw_data.min_mz * std::exp(fwhm_ref / mz_ref * i / mesh.k);
            } break;
            case Instrument::QUAD: {
                double delta_mz = (raw_data.max_mz - raw_data.min_mz) /
                                  static_cast<double>(mesh.n - 1);
                mesh.bins_mz[i] = raw_data.min_mz + delta_mz * i / mesh.k;
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
    {
        auto weights = std::vector<double>(n * m);
        for (size_t i = 0; i < raw_data.scans.size(); ++i) {
            const auto &scan = raw_data.scans[i];
            double current_rt = scan.retention_time;

            // Find the bin for the current retention time.
            size_t index_rt = y_index(raw_data, current_rt, mesh.t);

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
                size_t index_mz = x_index(raw_data, current_mz, mesh.k);

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

struct Peak {
    // ID of this peak. Should be kept for futher processing.
    size_t id;
    // Center of the peak in index space (Coordinates of local maxima).
    uint64_t local_max_i;
    uint64_t local_max_j;
    // Real mz/rt values for the center of this peak (From the local maxima
    // coordinates).
    double local_max_mz;
    double local_max_rt;
    // Height of the peak (Height of local maxima).
    double local_max_height;

    // Simple estimation of the peak metrics on the mesh values based on the
    // slope descent.
    //
    // Sumation of all intensities within the peak boundary. (Ignores holes,
    // i.e. does not interpolate values in case of non closed set).
    double slope_descent_total_intensity;
    // Estimated values for the position of the 2D peak based on the slope
    // descent points.
    double slope_descent_mz;
    double slope_descent_rt;
    // Estimated mz/rt values for the standard deviation of the peak in both
    // axes. (Ignores holes).
    double slope_descent_sigma_mz;
    double slope_descent_sigma_rt;
    // Average intensity on the boundary of the peak.
    double slope_descent_border_background;
    // NOTE: number of points within the boundary found via slope descent?

    // Region of interest for this peak.
    double roi_min_mz;
    double roi_max_mz;
    double roi_min_rt;
    double roi_max_rt;
    // Simple estimation of the peak metrics on the mesh values.
    double mesh_roi_mz;
    double mesh_roi_rt;
    double mesh_roi_sigma_mz;
    double mesh_roi_sigma_rt;
    double mesh_roi_total_intensity;
    // Simple estimation of the peak metrics on the raw data.
    double raw_roi_mz;
    double raw_roi_rt;
    double raw_roi_sigma_mz;
    double raw_roi_sigma_rt;
    double raw_roi_max_height;
    double raw_roi_total_intensity;
    double raw_roi_mean_height;
    double raw_roi_sigma_height;
    uint64_t raw_roi_num_points;
    uint64_t raw_roi_num_scans;

    // Matrix A for 2D gaussian fitting using least squares.
    double A[5][5];
    // Vector C for 2D gaussian fitting using least squares.
    double C[5];

    // FIXME: OMG WHAT IS HAPPENINGGGGGGGG
    double a_0_0() { return A[0][0]; }
    double a_0_1() { return A[0][1]; }
    double a_0_2() { return A[0][2]; }
    double a_0_3() { return A[0][3]; }
    double a_0_4() { return A[0][4]; }
    double a_1_0() { return A[1][0]; }
    double a_1_1() { return A[1][1]; }
    double a_1_2() { return A[1][2]; }
    double a_1_3() { return A[1][3]; }
    double a_1_4() { return A[1][4]; }
    double a_2_0() { return A[2][0]; }
    double a_2_1() { return A[2][1]; }
    double a_2_2() { return A[2][2]; }
    double a_2_3() { return A[2][3]; }
    double a_2_4() { return A[2][4]; }
    double a_3_0() { return A[3][0]; }
    double a_3_1() { return A[3][1]; }
    double a_3_2() { return A[3][2]; }
    double a_3_3() { return A[3][3]; }
    double a_3_4() { return A[3][4]; }
    double a_4_0() { return A[4][0]; }
    double a_4_1() { return A[4][1]; }
    double a_4_2() { return A[4][2]; }
    double a_4_3() { return A[4][3]; }
    double a_4_4() { return A[4][4]; }
    double c_0() { return C[0]; }
    double c_1() { return C[1]; }
    double c_2() { return C[2]; }
    double c_3() { return C[3]; }
    double c_4() { return C[4]; }
};

struct MeshIndex {
    size_t i;
    size_t j;
};

// FIXME: Probably this should be the default.
std::vector<MeshIndex> find_local_max_idx(const Mesh &mesh) {
    std::vector<MeshIndex> points;
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
                points.push_back({i, j});
            }
        }
    }

    auto sort_by_value = [&mesh](const MeshIndex &p1,
                                 const MeshIndex &p2) -> bool {
        return (mesh.matrix[p1.i + p1.j * mesh.n] >
                mesh.matrix[p2.i + p2.j * mesh.n]);
    };
    std::stable_sort(points.begin(), points.end(), sort_by_value);

    return points;
}

void explore_peak_slope(uint64_t i, uint64_t j, double previous_value,
                        const Mesh &mesh, std::vector<MeshIndex> &points) {
    // Check that the point has not being already included.
    for (const auto &point : points) {
        if (point.i == i && point.j == j) {
            return;
        }
    }

    double value = mesh.matrix[i + j * mesh.n];
    if (previous_value >= 0 && (previous_value < value || value <= 0.00001)) {
        return;
    }

    points.push_back({i, j});

    // Return if we are at the edge of the grid.
    if (i < 1 || i >= mesh.n - 1 || j < 1 || j >= mesh.m - 1) {
        return;
    }

    explore_peak_slope(i - 1, j, value, mesh, points);
    explore_peak_slope(i + 1, j, value, mesh, points);
    explore_peak_slope(i, j + 1, value, mesh, points);
    explore_peak_slope(i, j - 1, value, mesh, points);
    explore_peak_slope(i - 1, j - 1, value, mesh, points);
    explore_peak_slope(i + 1, j + 1, value, mesh, points);
    explore_peak_slope(i - 1, j + 1, value, mesh, points);
    explore_peak_slope(i + 1, j - 1, value, mesh, points);
}

std::vector<MeshIndex> find_boundary(std::vector<MeshIndex> &points) {
    // Under the constraints of the grid coordinates, we need at least 5 points
    // in order to have a boundary that does not contain all the points in the
    // initial set.
    if (points.size() < 5) {
        return points;
    }

    // Check if this point is a boundary by trying to find all 8 neighbours, if
    // the point does not have all of them, then it is a boundary point.
    auto point_exists = [&points](const MeshIndex &p) {
        for (const auto &point : points) {
            if (p.i == point.i && p.j == point.j) {
                return true;
            }
        }
        return false;
    };
    std::vector<MeshIndex> boundary;
    for (const auto &point : points) {
        if (!point_exists({point.i - 1, point.j - 1}) ||
            !point_exists({point.i, point.j - 1}) ||
            !point_exists({point.i + 1, point.j - 1}) ||
            !point_exists({point.i - 1, point.j}) ||
            !point_exists({point.i + 1, point.j}) ||
            !point_exists({point.i - 1, point.j + 1}) ||
            !point_exists({point.i, point.j + 1}) ||
            !point_exists({point.i + 1, point.j + 1})) {
            boundary.push_back(point);
        }
    }
    return boundary;
}

Peak build_peak(const RawData::RawData &raw_data, const Mesh &mesh,
                const MeshIndex &local_max) {
    Peak peak = {};
    peak.id = 0;
    peak.local_max_i = local_max.i;
    peak.local_max_j = local_max.j;
    peak.local_max_mz = mesh.bins_mz[local_max.i];
    peak.local_max_rt = mesh.bins_rt[local_max.j];
    peak.local_max_height = mesh.matrix[local_max.i + local_max.j * mesh.n];

    // Find the points within the boundary by slope descent on the mesh from the
    // local max.
    std::vector<MeshIndex> peak_points;
    // std::cout << peak.id << std::endl;
    explore_peak_slope(local_max.i, local_max.j, -1, mesh, peak_points);
    // FIXME: Should this just set NaN to boundary related peaks?
    if (peak_points.size() <= 1) {
        std::ostringstream error_stream;
        error_stream
            << "couldn't find any points on the mesh for this local max";
        throw std::invalid_argument(error_stream.str());
    }

    {
        // TODO(alex): error handling. What happens if the number of points is
        // very small? We should probably ignore peaks with less than 5 points
        // so that it has dimensionality in both mz and rt:
        //
        //   | |+| |
        //   |+|c|+|
        //   | |+| |
        std::vector<MeshIndex> peak_boundary;
        peak_boundary = find_boundary(peak_points);
        // FIXME: Should this just set NaN to boundary related peaks?
        if (peak_boundary.empty()) {
            std::ostringstream error_stream;
            error_stream << "couldn't find any points inside the boundary";
            throw std::invalid_argument(error_stream.str());
        }

        // Calculate the average background intensity from the boundary.
        double boundary_sum = 0;
        for (const auto &point : peak_boundary) {
            boundary_sum += mesh.matrix[point.i + point.j * mesh.n];
        }
        peak.slope_descent_border_background =
            boundary_sum / peak_boundary.size();
    }

    // Calculate the total ion intensity on the peak for the values on the grid
    // and the sigma in mz and rt. The sigma is calculated by using the
    // algebraic formula for the variance of the random variable X:
    //
    //     Var(X) = E[X^2] - E[X]
    //
    // Where E[X] is the estimated value for X.
    //
    // In order to generalize this formula for the 2D blob, all values at the
    // same index will be aggregated together.
    //
    // TODO(alex): Note that this can cause catastrophic cancellation or
    // loss of significance. Probably the best option is to use a variant of
    // the Welford's method for computing the variance in a single pass. See:
    //
    //     http://jonisalonen.com/2013/deriving-welfords-method-for-computing-variance/
    //     https://ipfs.io/ipfs/QmXoypizjW3WknFiJnKLwHCnL72vedxjQkDDP1mXWo6uco/wiki/Algorithms_for_calculating_variance.html
    //
    {
        double height_sum = 0;
        double x_sum = 0;
        double y_sum = 0;
        double x_sig = 0;
        double y_sig = 0;
        for (const auto &point : peak_points) {
            double mz = mesh.bins_mz[point.i];
            double rt = mesh.bins_rt[point.j];
            double value = mesh.matrix[point.i + point.j * mesh.n];

            height_sum += value;
            x_sum += value * mz;
            y_sum += value * rt;
            x_sig += value * mz * mz;
            y_sig += value * rt * rt;
        }
        peak.slope_descent_mz = x_sum / height_sum;
        peak.slope_descent_rt = y_sum / height_sum;
        peak.slope_descent_sigma_mz =
            std::sqrt((x_sig / height_sum) - std::pow(x_sum / height_sum, 2));
        peak.slope_descent_sigma_rt =
            std::sqrt((y_sig / height_sum) - std::pow(y_sum / height_sum, 2));
        peak.slope_descent_total_intensity = height_sum;
    }

    // Calculate the ROI for a given local max.
    {
        double mz = peak.local_max_mz;
        double rt = peak.local_max_rt;

        double theoretical_sigma_mz =
            fwhm_to_sigma(fwhm_at(raw_data, mesh.bins_mz[local_max.i]));
        double theoretical_sigma_rt = fwhm_to_sigma(raw_data.fwhm_rt);

        peak.roi_min_mz = mz - 3 * theoretical_sigma_mz;
        peak.roi_max_mz = mz + 3 * theoretical_sigma_mz;
        peak.roi_min_rt = rt - 3 * theoretical_sigma_rt;
        peak.roi_max_rt = rt + 3 * theoretical_sigma_rt;
    }

    // Calculate the estimation of values for the mesh points in the ROI.
    {
        // Find min_i via binary search.
        size_t min_i = 0;
        {
            size_t l = 0;
            size_t r = mesh.n - 1;
            while (l <= r) {
                min_i = (l + r) / 2;
                if (mesh.bins_mz[min_i] < peak.roi_min_mz) {
                    l = min_i + 1;
                } else if (mesh.bins_mz[min_i] > peak.roi_min_mz) {
                    r = min_i - 1;
                } else {
                    break;
                }
                if (min_i == 0) {
                    break;
                }
            }
        }
        // Find min_j via binary search.
        size_t min_j = 0;
        {
            size_t l = 0;
            size_t r = mesh.m - 1;
            while (l <= r) {
                min_j = (l + r) / 2;
                if (mesh.bins_rt[min_j] < peak.roi_min_rt) {
                    l = min_j + 1;
                } else if (mesh.bins_rt[min_j] > peak.roi_min_rt) {
                    r = min_j - 1;
                } else {
                    break;
                }
                if (min_j == 0) {
                    break;
                }
            }
        }
        double height_sum = 0;
        double x_sum = 0;
        double y_sum = 0;
        double x_sig = 0;
        double y_sig = 0;
        for (size_t j = min_j; j < mesh.m; ++j) {
            if (mesh.bins_rt[j] > peak.roi_max_rt) {
                break;
            }
            for (size_t i = min_i; i < mesh.n; ++i) {
                if (mesh.bins_mz[i] > peak.roi_max_mz) {
                    break;
                }
                double mz = mesh.bins_mz[i];
                double rt = mesh.bins_rt[j];
                double value = mesh.matrix[i + j * mesh.n];
                height_sum += value;
                x_sum += value * mz;
                y_sum += value * rt;
                x_sig += value * mz * mz;
                y_sig += value * rt * rt;
            }
        }
        peak.mesh_roi_mz = x_sum / height_sum;
        peak.mesh_roi_rt = y_sum / height_sum;
        peak.mesh_roi_sigma_mz =
            std::sqrt((x_sig / height_sum) - std::pow(x_sum / height_sum, 2));
        peak.mesh_roi_sigma_rt =
            std::sqrt((y_sig / height_sum) - std::pow(y_sum / height_sum, 2));
        peak.mesh_roi_total_intensity = height_sum;
    }

    {
        // Initializing A and C
        // Matrix A.
        // Row 0
        peak.A[0][0] = 0;
        peak.A[0][1] = 0;
        peak.A[0][2] = 0;
        peak.A[0][3] = 0;
        peak.A[0][4] = 0;
        // Row 1
        peak.A[1][0] = 0;
        peak.A[1][1] = 0;
        peak.A[1][2] = 0;
        peak.A[1][3] = 0;
        peak.A[1][4] = 0;
        // Row 2
        peak.A[2][0] = 0;
        peak.A[2][1] = 0;
        peak.A[2][2] = 0;
        peak.A[2][3] = 0;
        peak.A[2][4] = 0;
        // Row 3
        peak.A[3][0] = 0;
        peak.A[3][1] = 0;
        peak.A[3][2] = 0;
        peak.A[3][3] = 0;
        peak.A[3][4] = 0;
        // Row 4
        peak.A[4][0] = 0;
        peak.A[4][1] = 0;
        peak.A[4][2] = 0;
        peak.A[4][3] = 0;
        peak.A[4][4] = 0;
        // Vector C.
        peak.C[0] = 0;
        peak.C[1] = 0;
        peak.C[2] = 0;
        peak.C[3] = 0;
        peak.C[4] = 0;

        const auto &scans = raw_data.scans;
        // FIXME: Make nan instead?
        if (scans.size() == 0) {
            std::ostringstream error_stream;
            error_stream << "the given raw_data is empty";
            throw std::invalid_argument(error_stream.str());
        }

        // Find scan indices.
        double internal_min_rt = peak.roi_min_rt;
        double internal_max_rt = peak.roi_max_rt;
        if (internal_min_rt < raw_data.min_rt) {
            internal_min_rt = raw_data.min_rt;
        }
        if (internal_max_rt > raw_data.max_rt) {
            internal_max_rt = raw_data.max_rt;
        }

        // Binary search for lower rt bound.
        size_t min_j = 0;
        size_t max_j = scans.size();
        size_t l = min_j;
        size_t r = max_j - 1;
        while (l <= r) {
            min_j = (l + r) / 2;
            if (scans[min_j].retention_time < internal_min_rt) {
                l = min_j + 1;
            } else if (scans[min_j].retention_time > internal_min_rt) {
                r = min_j - 1;
            } else {
                break;
            }
            if (min_j == 0) {
                break;
            }
        }
        double height_sum = 0;
        double sq_height_sum = 0;
        double x_sum = 0;
        double y_sum = 0;
        double x_sig = 0;
        double y_sig = 0;
        for (size_t j = min_j; j < max_j; ++j) {
            const auto &scan = scans[j];
            if (scan.num_points == 0) {
                continue;
            }
            if (scan.retention_time > internal_max_rt) {
                break;
            }
            ++peak.raw_roi_num_scans;
            // Binary search for lower mz bound.
            double internal_min_mz = peak.roi_min_mz;
            double internal_max_mz = peak.roi_max_mz;
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
                double mz = scan.mz[i];
                double rt = scan.retention_time;
                double value = scan.intensity[i];
                if (value > peak.raw_roi_max_height) {
                    peak.raw_roi_max_height = value;
                }
                ++peak.raw_roi_num_points;
                height_sum += value;
                sq_height_sum += value * value;
                x_sum += value * mz;
                y_sum += value * rt;
                x_sig += value * mz * mz;
                y_sig += value * rt * rt;
            }
        }
        // FIXME: Not controlling for div/0.
        peak.raw_roi_mz = x_sum / height_sum;
        peak.raw_roi_rt = y_sum / height_sum;
        peak.raw_roi_sigma_mz =
            std::sqrt((x_sig / height_sum) - std::pow(x_sum / height_sum, 2));
        peak.raw_roi_sigma_rt =
            std::sqrt((y_sig / height_sum) - std::pow(y_sum / height_sum, 2));
        peak.raw_roi_total_intensity = height_sum;
        peak.raw_roi_mean_height = height_sum / peak.raw_roi_num_points;
        peak.raw_roi_sigma_height =
            std::sqrt(sq_height_sum / peak.raw_roi_num_points -
                      peak.raw_roi_mean_height * peak.raw_roi_mean_height);

        for (size_t j = min_j; j < max_j; ++j) {
            const auto &scan = scans[j];
            if (scan.num_points == 0) {
                continue;
            }
            if (scan.retention_time > internal_max_rt) {
                break;
            }
            // Binary search for lower mz bound.
            double internal_min_mz = peak.roi_min_mz;
            double internal_max_mz = peak.roi_max_mz;
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
                double mz = scan.mz[i];
                double rt = scan.retention_time;
                double value = scan.intensity[i];

                // Calculate the values for the A matrix and C vector necessary
                // for the 2D Gaussian fitting using least squares.
                {
                    // FIXME: We might need to center the mz/rt values.
                    double x = (mz - peak.raw_roi_mz);
                    double y = (rt - peak.raw_roi_rt);
                    double z = value;
                    double log_z = std::log(z);
                    double z_2 = std::pow(z, 2);
                    double x_2 = std::pow(x, 2);
                    double x_3 = std::pow(x, 3);
                    double x_4 = std::pow(x, 4);
                    double y_2 = std::pow(y, 2);
                    double y_3 = std::pow(y, 3);
                    double y_4 = std::pow(y, 4);
                    // Matrix A.
                    // Row 0
                    peak.A[0][0] += z_2;
                    peak.A[0][1] += z_2 * x;
                    peak.A[0][2] += z_2 * x_2;
                    peak.A[0][3] += z_2 * y;
                    peak.A[0][4] += z_2 * y_2;
                    // Row 1
                    peak.A[1][0] += z_2 * x;
                    peak.A[1][1] += z_2 * x_2;
                    peak.A[1][2] += z_2 * x_3;
                    peak.A[1][3] += z_2 * x * y;
                    peak.A[1][4] += z_2 * x * y_2;
                    // Row 2
                    peak.A[2][0] += z_2 * x_2;
                    peak.A[2][1] += z_2 * x_3;
                    peak.A[2][2] += z_2 * x_4;
                    peak.A[2][3] += z_2 * x_2 * y;
                    peak.A[2][4] += z_2 * x_2 * y_2;
                    // Row 3
                    peak.A[3][0] += z_2 * y;
                    peak.A[3][1] += z_2 * x * y;
                    peak.A[3][2] += z_2 * x_2 * y;
                    peak.A[3][3] += z_2 * y_2;
                    peak.A[3][4] += z_2 * y_3;
                    // Row 4
                    peak.A[4][0] += z_2 * y_2;
                    peak.A[4][1] += z_2 * x * y_2;
                    peak.A[4][2] += z_2 * x_2 * y_2;
                    peak.A[4][3] += z_2 * y_3;
                    peak.A[4][4] += z_2 * y_4;
                    // Vector C.
                    peak.C[0] += z_2 * log_z;
                    peak.C[1] += z_2 * x * log_z;
                    peak.C[2] += z_2 * x_2 * log_z;
                    peak.C[3] += z_2 * y * log_z;
                    peak.C[4] += z_2 * y_2 * log_z;
                }
            }
        }

        // FIXME: Make nan instead?
        // if (raw_points.num_points == 0) {
        // std::ostringstream error_stream;
        // error_stream << "couldn't find raw_data points on the given ROI";
        // throw std::invalid_argument(error_stream.str());
        //}
    }

    return peak;
}

std::vector<Peak> find_peaks(const RawData::RawData &raw_data,
                             const Mesh &mesh) {
    auto local_max = find_local_max_idx(mesh);
    std::vector<Peak> peaks;
    for (const auto &lm : local_max) {
        peaks.push_back(build_peak(raw_data, mesh, lm));
    }
    // TODO: Sort peaks
    // TODO: Fill peak ids.
    return peaks;
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

std::vector<std::vector<Peak>> warp_peaks(
    std::vector<std::vector<Peak>> all_peaks, size_t reference_index,
    int64_t slack, int64_t window_size, int64_t num_points,
    double rt_expand_factor, int64_t peaks_per_window) {
    // TODO(alex): Validate the parameters and throw an error if appropriate.
    Warp2D::Parameters parameters = {slack, window_size, num_points,
                                     peaks_per_window, rt_expand_factor};
    std::vector<std::vector<Peak>> all_warped_peaks;
    const auto &reference_peaks_original = all_warped_peaks[reference_index];
    auto translate_peak_format = [](const std::vector<Peak> &before_peaks)
        -> std::vector<Centroid::Peak> {
        auto after_peaks = std::vector<Centroid::Peak>(before_peaks.size());
        for (size_t i = 0; i < before_peaks.size(); ++i) {
            after_peaks[i].i = before_peaks[i].local_max_i;
            after_peaks[i].j = before_peaks[i].local_max_j;
            after_peaks[i].mz = before_peaks[i].local_max_mz;
            after_peaks[i].rt = before_peaks[i].local_max_rt;
            after_peaks[i].height = before_peaks[i].local_max_height;

            // NOTE: Currently using slope_descent quantification.
            after_peaks[i].total_intensity =
                before_peaks[i].slope_descent_total_intensity;
            after_peaks[i].sigma_mz = before_peaks[i].slope_descent_sigma_mz;
            after_peaks[i].sigma_rt = before_peaks[i].slope_descent_sigma_rt;
            after_peaks[i].border_background =
                before_peaks[i].slope_descent_border_background;
        }
        return after_peaks;
    };
    auto reference_peaks = translate_peak_format(reference_peaks_original);

    for (size_t i = 0; i < all_peaks.size(); ++i) {
        if (i == reference_index) {
            continue;
        }
        auto peaks = translate_peak_format(all_peaks[i]);
        std::vector<Centroid::Peak> warped_peaks;
        warped_peaks =
            Warp2D::Runners::Serial::run(reference_peaks, peaks, parameters);
    }
    return all_warped_peaks;
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

    py::class_<PythonAPI::Peak>(m, "Peak")
        .def_readonly("id", &PythonAPI::Peak::id)
        .def_readonly("local_max_i", &PythonAPI::Peak::local_max_i)
        .def_readonly("local_max_j", &PythonAPI::Peak::local_max_j)
        .def_readonly("local_max_mz", &PythonAPI::Peak::local_max_mz)
        .def_readonly("local_max_rt", &PythonAPI::Peak::local_max_rt)
        .def_readonly("local_max_height", &PythonAPI::Peak::local_max_height)
        .def_readonly("slope_descent_mz", &PythonAPI::Peak::slope_descent_mz)
        .def_readonly("slope_descent_rt", &PythonAPI::Peak::slope_descent_rt)
        .def_readonly("slope_descent_sigma_mz",
                      &PythonAPI::Peak::slope_descent_sigma_mz)
        .def_readonly("slope_descent_sigma_rt",
                      &PythonAPI::Peak::slope_descent_sigma_rt)
        .def_readonly("slope_descent_total_intensity",
                      &PythonAPI::Peak::slope_descent_total_intensity)
        .def_readonly("slope_descent_border_background",
                      &PythonAPI::Peak::slope_descent_border_background)
        .def_readonly("roi_min_mz", &PythonAPI::Peak::roi_min_mz)
        .def_readonly("roi_max_mz", &PythonAPI::Peak::roi_max_mz)
        .def_readonly("roi_min_rt", &PythonAPI::Peak::roi_min_rt)
        .def_readonly("roi_max_rt", &PythonAPI::Peak::roi_max_rt)
        .def_readonly("mesh_roi_mz", &PythonAPI::Peak::mesh_roi_mz)
        .def_readonly("mesh_roi_rt", &PythonAPI::Peak::mesh_roi_rt)
        .def_readonly("mesh_roi_sigma_mz", &PythonAPI::Peak::mesh_roi_sigma_mz)
        .def_readonly("mesh_roi_sigma_rt", &PythonAPI::Peak::mesh_roi_sigma_rt)
        .def_readonly("mesh_roi_total_intensity",
                      &PythonAPI::Peak::mesh_roi_total_intensity)
        .def_readonly("raw_roi_mz", &PythonAPI::Peak::raw_roi_mz)
        .def_readonly("raw_roi_rt", &PythonAPI::Peak::raw_roi_rt)
        .def_readonly("raw_roi_sigma_mz", &PythonAPI::Peak::raw_roi_sigma_mz)
        .def_readonly("raw_roi_sigma_rt", &PythonAPI::Peak::raw_roi_sigma_rt)
        .def_readonly("raw_roi_total_intensity",
                      &PythonAPI::Peak::raw_roi_total_intensity)
        .def_readonly("raw_roi_max_height",
                      &PythonAPI::Peak::raw_roi_max_height)
        .def_readonly("raw_roi_mean_height",
                      &PythonAPI::Peak::raw_roi_mean_height)
        .def_readonly("raw_roi_sigma_height",
                      &PythonAPI::Peak::raw_roi_sigma_height)
        .def_readonly("raw_roi_num_points",
                      &PythonAPI::Peak::raw_roi_num_points)
        .def_readonly("raw_roi_num_scans", &PythonAPI::Peak::raw_roi_num_scans)
        .def("a_0_0", &PythonAPI::Peak::a_0_0)
        .def("a_0_1", &PythonAPI::Peak::a_0_1)
        .def("a_0_2", &PythonAPI::Peak::a_0_2)
        .def("a_0_3", &PythonAPI::Peak::a_0_3)
        .def("a_0_4", &PythonAPI::Peak::a_0_4)
        .def("a_1_0", &PythonAPI::Peak::a_1_0)
        .def("a_1_1", &PythonAPI::Peak::a_1_1)
        .def("a_1_2", &PythonAPI::Peak::a_1_2)
        .def("a_1_3", &PythonAPI::Peak::a_1_3)
        .def("a_1_4", &PythonAPI::Peak::a_1_4)
        .def("a_2_0", &PythonAPI::Peak::a_2_0)
        .def("a_2_1", &PythonAPI::Peak::a_2_1)
        .def("a_2_2", &PythonAPI::Peak::a_2_2)
        .def("a_2_3", &PythonAPI::Peak::a_2_3)
        .def("a_2_4", &PythonAPI::Peak::a_2_4)
        .def("a_3_0", &PythonAPI::Peak::a_3_0)
        .def("a_3_1", &PythonAPI::Peak::a_3_1)
        .def("a_3_2", &PythonAPI::Peak::a_3_2)
        .def("a_3_3", &PythonAPI::Peak::a_3_3)
        .def("a_3_4", &PythonAPI::Peak::a_3_4)
        .def("a_4_0", &PythonAPI::Peak::a_4_0)
        .def("a_4_1", &PythonAPI::Peak::a_4_1)
        .def("a_4_2", &PythonAPI::Peak::a_4_2)
        .def("a_4_3", &PythonAPI::Peak::a_4_3)
        .def("a_4_4", &PythonAPI::Peak::a_4_4)
        .def("c_0", &PythonAPI::Peak::c_0)
        .def("c_1", &PythonAPI::Peak::c_1)
        .def("c_2", &PythonAPI::Peak::c_2)
        .def("c_3", &PythonAPI::Peak::c_3)
        .def("c_4", &PythonAPI::Peak::c_4);

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
             py::arg("max_rt"))
        .def("find_peaks", &PythonAPI::find_peaks,
             "Find all peaks in the given mesh", py::arg("raw_data"),
             py::arg("mesh"))
        .def("warp_peaks", &PythonAPI::warp_peaks,
             "Warp peak lists to maximize the similarity with the given "
             "reference",
             py::arg("all_peaks"), py::arg("reference_index"), py::arg("slack"),
             py::arg("window_size"), py::arg("num_points"),
             py::arg("rt_expand_factor"), py::arg("peaks_per_window"));
}
