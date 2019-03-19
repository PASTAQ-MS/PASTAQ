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

// TODO: This should be tested.
RoiIndices find_roi_indexes(const RawData::RawData &raw_data, const Roi &roi) {
    // The search procedure for the indices for both retention time and m/z is
    // the same: Fist a binary search is performed to find the lower bound,
    // followed by a linear search from the minimum index to find the upper
    // bound.
    //
    const auto &scans = raw_data.scans;
    if (scans.size() == 0) {
        // return {0, 0};
        // TODO: Throw exception for python code to catch?
    }
    RoiIndices roi_indices;
    // Find scan indices.
    {
        double min_rt = roi.min_rt;
        double max_rt = roi.max_rt;
        if (min_rt < raw_data.min_rt) {
            min_rt = raw_data.min_rt;
        }
        if (max_rt > raw_data.max_rt) {
            max_rt = raw_data.max_rt;
        }
        // Binary search for lower bound.
        size_t min_j = 0;
        size_t max_j = raw_data.scans.size() - 1;
        size_t l = min_j;
        size_t r = max_j;
        while (l <= r) {
            min_j = (l + r) / 2;
            if (scans[min_j].retention_time < min_rt) {
                l = min_j + 1;
            } else if (scans[min_j].retention_time > min_rt) {
                r = min_j - 1;
            } else {
                break;
            }
        }
        // Linear search for upper bound.
        if (min_j != max_j) {
            min_j = min_j + 1;
        }
        for (size_t i = min_j; i < scans.size(); ++i) {
            if (scans[i].retention_time >= max_rt) {
                max_j = i;
                break;
            }
        }
        roi_indices.min_j = min_j;
        roi_indices.max_j = max_j;
    }

    // Find m/z indices per scan.
    roi_indices.min_i =
        std::vector<size_t>(roi_indices.max_j - roi_indices.min_j + 1);
    roi_indices.max_i =
        std::vector<size_t>(roi_indices.max_j - roi_indices.min_j + 1);
    size_t k = 0;
    for (size_t j = roi_indices.min_j; j < roi_indices.max_j; ++j, ++k) {
        const auto &scan = scans[j];
        if (scan.num_points == 0) {
            // return {0, 0};
            // TODO: Throw exception for python code to catch?
            continue;
        }
        {
            double min_mz = roi.min_mz;
            double max_mz = roi.max_mz;
            if (min_mz < scan.mz[0]) {
                min_mz = scan.mz[0];
            }
            if (max_mz > scan.mz[scan.num_points - 1]) {
                max_mz = scan.mz[scan.num_points - 1];
            }
            // Binary search for lower bound.
            size_t min_i = 0;
            size_t max_i = scan.num_points - 1;
            size_t l = min_i;
            size_t r = max_i;
            while (l <= r) {
                min_i = (l + r) / 2;
                if (scan.mz[min_i] < min_mz) {
                    l = min_i + 1;
                } else if (scan.mz[min_i] > min_mz) {
                    r = min_i - 1;
                } else {
                    break;
                }
            }
            // Linear search for upper bound.
            if (min_i != max_i) {
                min_i = min_i + 1;
            }
            for (size_t i = min_i; i < scan.mz.size(); ++i) {
                if (scan.mz[i] >= max_mz) {
                    max_i = i;
                    break;
                }
            }
            roi_indices.min_i[k] = min_i;
            roi_indices.max_i[k] = max_i;
        }
    }

    return roi_indices;
}

struct RawPoint {
    double rt;
    double mz;
    double intensity;
};

struct RawPoints {
    std::vector<double> rt;
    std::vector<double> mz;
    std::vector<double> intensity;
};

RawPoints find_raw_points(const RawData::RawData &raw_data, double min_mz,
                          double max_mz, double min_rt, double max_rt) {
    auto indices = find_roi_indexes(raw_data, {min_mz, max_mz, min_rt, max_rt});
    RawPoints raw_points;
    size_t k = 0;
    if (((int64_t)indices.max_j - (int64_t)indices.min_j) <= 3) {
        std::ostringstream error_stream;
        error_stream << "not enough scans to perform fitting...";
        throw std::invalid_argument(error_stream.str());
    }
    for (size_t j = indices.min_j; j < indices.max_j; ++j, ++k) {
        const auto &scan = raw_data.scans[j];
        if (((int64_t)indices.max_i[k] - (int64_t)indices.min_i[k]) <= 3) {
            continue;
        }
        for (size_t i = indices.min_i[k]; i < indices.max_i[k]; ++i) {
            raw_points.rt.push_back(scan.retention_time);
            raw_points.mz.push_back(scan.mz[i]);
            raw_points.intensity.push_back(scan.intensity[i]);
        }
    }
    return raw_points;
}

uint64_t x_index(const RawData::RawData &raw_data, double mz, uint64_t k) {
    // FIXME: This only works for ORBITRAP data for now.
    double fwhm_ref = raw_data.reference_mz / raw_data.resolution_ms1;
    double a = fwhm_ref / std::pow(raw_data.reference_mz, 1.5);
    double b = (1 / std::sqrt(raw_data.min_mz) - 1 / std::sqrt(mz));
    return k * 2 / a * b;
}

uint64_t y_index(const RawData::RawData &raw_data, double rt, uint64_t k) {
    double delta_rt = raw_data.fwhm_rt / k;
    return std::ceil((rt - raw_data.min_rt) / delta_rt);
}

std::tuple<uint64_t, uint64_t> calculate_dimensions(
    const RawData::RawData &raw_data, uint64_t num_samples_per_peak_mz,
    uint64_t num_samples_per_peak_rt) {
    // FIXME: This only works for ORBITRAP data for now.
    return std::tuple<uint64_t, uint64_t>(
        x_index(raw_data, raw_data.max_mz, num_samples_per_peak_mz) + 1,
        y_index(raw_data, raw_data.max_rt, num_samples_per_peak_rt) + 1);
}

double mz_at(const RawData::RawData &raw_data, uint64_t num_samples_per_peak_mz,
             uint64_t n) {
    // FIXME: This only works for ORBITRAP data for now.
    double a = 1 / std::sqrt(raw_data.min_mz);
    double fwhm_ref = raw_data.reference_mz / raw_data.resolution_ms1;
    double b = fwhm_ref / std::pow(raw_data.reference_mz, 1.5) * n / 2 /
               num_samples_per_peak_mz;
    double c = a - b;
    return 1 / (c * c);
}

double fwhm_at(const RawData::RawData &raw_data, double mz) {
    // FIXME: This only works for ORBITRAP data for now.
    double fwhm_ref = raw_data.reference_mz / raw_data.resolution_ms1;
    return fwhm_ref * std::pow(mz / raw_data.reference_mz, 1.5);
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
        // FIXME: For now...
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

// The given smoothing coefficients will determine how much smoothing we apply
// on a given mesh. Since the smoothing is calculated as a factor of a Gaussian
// sigma, the smoothing ratio will be the same regardless of the sampling size.
// Smoothing ratios of 0.5-1 are recommended as a starting point.
Mesh resample(const RawData::RawData &raw_data,
              uint64_t num_samples_per_peak_mz,
              uint64_t num_samples_per_peak_rt, double smoothing_coef_mz,
              double smoothing_coef_rt) {
    auto [n, m] = calculate_dimensions(raw_data, num_samples_per_peak_mz,
                                       num_samples_per_peak_rt);
    Mesh mesh;
    mesh.n = n;
    mesh.m = m;
    mesh.matrix = std::vector<double>(n * m);
    mesh.bins_mz = std::vector<double>(n);
    mesh.bins_rt = std::vector<double>(m);

    // Generate bins_mz.
    for (size_t i = 0; i < n; ++i) {
        mesh.bins_mz[i] = mz_at(raw_data, num_samples_per_peak_mz, i);
    }
    // Generate bins_rt.
    double delta_rt = (raw_data.max_rt - raw_data.min_rt) / (m - 1);
    for (size_t j = 0; j < m; ++j) {
        mesh.bins_rt[j] = raw_data.min_rt + delta_rt * j;
    }

    // Pre-calculate the smoothing sigma values for all bins of the grid.
    double sigma_rt = fwhm_to_sigma(raw_data.fwhm_rt) * smoothing_coef_rt;
    auto sigma_mz_vect = std::vector<double>(n);
    for (size_t i = 0; i < n; ++i) {
        sigma_mz_vect[i] = fwhm_to_sigma(fwhm_at(raw_data, mesh.bins_mz[i])) *
                           smoothing_coef_mz;
    }

    for (const auto &scan : raw_data.scans) {
        // Calculate the min and max indexes for retention time.

        // Find the bin for the current retention time.
        double current_rt = scan.retention_time;
        size_t index_rt =
            y_index(raw_data, current_rt, num_samples_per_peak_rt);

        // The smoothing kernel in rt is +-(num_samples_per_peak_rt/2).
        // int64_t j_min = index_rt - num_samples_per_peak_rt;
        int64_t j_min = y_index(raw_data, current_rt - 3 * sigma_rt,
                                num_samples_per_peak_rt);
        if (j_min < 0) {
            j_min = 0;
        }
        // int64_t j_max = index_rt + num_samples_per_peak_rt;
        int64_t j_max = y_index(raw_data, current_rt + 3 * sigma_rt,
                                num_samples_per_peak_rt);
        if (j_max >= m) {
            j_max = m - 1;
        }

        for (size_t k = 0; k < scan.num_points; ++k) {
            double current_intensity = scan.intensity[k];

            // Find the bin for the current mz.
            double current_mz = scan.mz[k];

            size_t index_mz =
                x_index(raw_data, current_mz, num_samples_per_peak_mz);
            double sigma_mz = sigma_mz_vect[index_mz];
            int64_t i_min = x_index(raw_data, current_mz - 3 * sigma_mz,
                                    num_samples_per_peak_mz);
            if (i_min < 0) {
                i_min = 0;
            }
            int64_t i_max = x_index(raw_data, current_mz + 3 * sigma_mz,
                                    num_samples_per_peak_mz);
            if (i_max >= n) {
                i_max = n - 1;
            }

            for (size_t j = j_min; j <= j_max; ++j) {
                for (size_t i = i_min; i <= i_max; ++i) {
                    // No need to do boundary check, since we are sure we are
                    // inside the grid.
                    double x = mesh.bins_mz[i];
                    double y = mesh.bins_rt[j];

                    // Calculate the gaussian weight for this point.
                    double a = (x - current_mz) / sigma_mz;
                    double b = (y - current_rt) / sigma_rt;
                    double weight = std::exp(-0.5 * (a * a + b * b));

                    // Set the value, weight and counts.
                    mesh.matrix[i + j * n] += weight * current_intensity;
                }
            }
        }
    }
    return mesh;
}

// The given smoothing coefficients will determine how much smoothing we apply
// on a given mesh. Since the smoothing is calculated as a factor of a Gaussian
// sigma, the smoothing ratio will be the same regardless of the sampling size.
// Smoothing ratios of 0.5-1 are recommended as a starting point.
//
// NOTE: Here we are trying an updated gaussian smoothing algorithm.
Mesh resample_2(const RawData::RawData &raw_data,
                uint64_t num_samples_per_peak_mz,
                uint64_t num_samples_per_peak_rt, double smoothing_coef_mz,
                double smoothing_coef_rt) {
    auto [n, m] = calculate_dimensions(raw_data, num_samples_per_peak_mz,
                                       num_samples_per_peak_rt);
    Mesh mesh;
    mesh.n = n;
    mesh.m = m;
    mesh.matrix = std::vector<double>(n * m);
    mesh.bins_mz = std::vector<double>(n);
    mesh.bins_rt = std::vector<double>(m);

    // Generate bins_mz.
    for (size_t i = 0; i < n; ++i) {
        mesh.bins_mz[i] = mz_at(raw_data, num_samples_per_peak_mz, i);
    }
    // Generate bins_rt.
    double delta_rt = (raw_data.max_rt - raw_data.min_rt) / (m - 1);
    for (size_t j = 0; j < m; ++j) {
        mesh.bins_rt[j] = raw_data.min_rt + delta_rt * j;
    }

    // Pre-calculate the smoothing sigma values for all bins of the grid.
    double sigma_rt = fwhm_to_sigma(raw_data.fwhm_rt) * smoothing_coef_rt;
    auto sigma_mz_vect = std::vector<double>(n);
    for (size_t i = 0; i < n; ++i) {
        sigma_mz_vect[i] = fwhm_to_sigma(fwhm_at(raw_data, mesh.bins_mz[i])) *
                           smoothing_coef_mz;
    }

    for (size_t row = 0; row < m - 1; ++row) {
        // Find scans that belong to this bin in rt.
        double min_rt = mesh.bins_rt[row] - 3 * sigma_rt;
        double max_rt = mesh.bins_rt[row] + 3 * sigma_rt;
        double current_rt = mesh.bins_rt[row];

        auto smoothed_row = std::vector<double>(n);
        for (size_t col = 0; col < n; ++col) {
            double min_mz = mesh.bins_mz[col] - 3 * sigma_mz_vect[col];
            double max_mz = mesh.bins_mz[col] + 3 * sigma_mz_vect[col];
            double current_mz = mesh.bins_mz[col];
            double sigma_mz = sigma_mz_vect[col];

            auto indices =
                find_roi_indexes(raw_data, {min_mz, max_mz, min_rt, max_rt});

            size_t k = 0;
            double sum_weights = 0;
            double sum_weights_values = 0;
            for (size_t j = indices.min_j; j < indices.max_j; ++j, ++k) {
                const auto &scan = raw_data.scans[j];
                if (((int64_t)indices.max_i[k] - (int64_t)indices.min_i[k]) ==
                    0) {
                    continue;
                }
                for (size_t i = indices.min_i[k]; i < indices.max_i[k]; ++i) {
                    double a = (scan.mz[i] - current_mz) / sigma_mz;
                    double b = (scan.retention_time - current_rt) / sigma_rt;
                    double weight = std::exp(-0.5 * (a * a + b * b));
                    sum_weights += weight;
                    sum_weights_values += weight * scan.intensity[i];
                }
            }

            // smoothed_row[col] = sum_weights_values / (sum_weights + 1);
            if (sum_weights == 0) {
                sum_weights = 1;
            }
            smoothed_row[col] = sum_weights_values / (sum_weights);
            mesh.matrix[col + row * n] = smoothed_row[col];
        }
    }

    return mesh;
}

// NOTE: Same as resample_3 but done in 2 passes.
Mesh resample_3(const RawData::RawData &raw_data,
                uint64_t num_samples_per_peak_mz,
                uint64_t num_samples_per_peak_rt, double smoothing_coef_mz,
                double smoothing_coef_rt) {
    auto [n, m] = calculate_dimensions(raw_data, num_samples_per_peak_mz,
                                       num_samples_per_peak_rt);
    Mesh mesh;
    mesh.n = n;
    mesh.m = m;
    mesh.matrix = std::vector<double>(n * m);
    mesh.bins_mz = std::vector<double>(n);
    mesh.bins_rt = std::vector<double>(m);

    // Generate bins_mz.
    for (size_t i = 0; i < n; ++i) {
        mesh.bins_mz[i] = mz_at(raw_data, num_samples_per_peak_mz, i);
    }
    // Generate bins_rt.
    double delta_rt = (raw_data.max_rt - raw_data.min_rt) / (m - 1);
    for (size_t j = 0; j < m; ++j) {
        mesh.bins_rt[j] = raw_data.min_rt + delta_rt * j;
    }

    // Pre-calculate the smoothing sigma values for all bins of the grid.
    double sigma_rt = fwhm_to_sigma(raw_data.fwhm_rt) * smoothing_coef_rt;
    auto sigma_mz_vect = std::vector<double>(n);
    for (size_t i = 0; i < n; ++i) {
        sigma_mz_vect[i] = fwhm_to_sigma(fwhm_at(raw_data, mesh.bins_mz[i])) *
                           smoothing_coef_mz;
    }

    std::vector<std::vector<double>> smoothed_scans;

    // MZ Smoothing.
    for (const auto &scan : raw_data.scans) {
        auto mzs = scan.mz;
        auto intensities = scan.intensity;
        size_t previous_min_i = 0;
        auto smoothed_row = std::vector<double>(n);
        for (size_t col = 0; col < n; ++col) {
            double current_mz = mesh.bins_mz[col];
            double sigma_mz = sigma_mz_vect[col];
            double min_mz = current_mz - 3 * sigma_mz;
            double max_mz = current_mz + 3 * sigma_mz;

            // Find mz points on the raw_data to be used for smoothing.
            int64_t min_i = 0;
            int64_t max_i = scan.num_points;
            for (size_t i = previous_min_i; i < scan.num_points; ++i) {
                if (scan.mz[i] > min_mz && i != 0) {
                    min_i = i - 1;
                    break;
                }
            }
            for (size_t i = min_i; i < scan.num_points; ++i) {
                if (scan.mz[i] > max_mz) {
                    max_i = i;
                    break;
                }
            }
            previous_min_i = min_i;

            // Get the smoothed mz estimate for this point.
            double sum_weights = 0;
            double sum_weights_values = 0;
            for (size_t i = min_i; i < max_i; ++i) {
                double a = (scan.mz[i] - current_mz) / sigma_mz;
                double weight = std::exp(-0.5 * (a * a));
                sum_weights += weight;
                sum_weights_values += weight * scan.intensity[i];
            }
            if (sum_weights == 0) {
                sum_weights = 1;
            }
            smoothed_row[col] = sum_weights_values / sum_weights;
        }
        smoothed_scans.emplace_back(smoothed_row);
    }

    // RT Smoothing.
    for (size_t row = 0; row < m; ++row) {
        double current_rt = mesh.bins_rt[row];
        double min_rt = current_rt - 3 * sigma_rt;
        double max_rt = current_rt + 3 * sigma_rt;

        // Find index range on the smoothed raw spectra for rt smoothing.
        size_t previous_min_j = 0;
        int64_t min_j = 0;
        int64_t max_j = raw_data.scans.size();
        for (size_t j = previous_min_j; j < raw_data.scans.size(); ++j) {
            if (raw_data.scans[j].retention_time > min_rt && j != 0) {
                min_j = j - 1;
                break;
            }
        }
        for (size_t j = min_j; j < raw_data.scans.size(); ++j) {
            if (raw_data.scans[j].retention_time > max_rt) {
                max_j = j;
                break;
            }
        }
        previous_min_j = min_j;

        for (size_t col = 0; col < n; ++col) {
            // Get the smoothed rt estimate for this point.
            double sum_weights = 0;
            double sum_weights_values = 0;
            for (size_t j = min_j; j < max_j; ++j) {
                double a =
                    (raw_data.scans[j].retention_time - current_rt) / sigma_rt;
                double weight = std::exp(-0.5 * (a * a));
                sum_weights += weight;
                sum_weights_values += weight * smoothed_scans[j][col];
            }
            if (sum_weights == 0) {
                sum_weights = 1;
            }

            mesh.matrix[col + row * n] = sum_weights_values / sum_weights;
        }
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
        .def("mz_at", &PythonAPI::mz_at,
             "Calculate the mz at the given N for the given raw file",
             py::arg("raw_data"), py::arg("num_mz") = 10, py::arg("n"))
        .def("fwhm_at", &PythonAPI::fwhm_at,
             "Calculate the width of the peak at the given m/z for the given "
             "raw file",
             py::arg("raw_data"), py::arg("mz"))
        .def("resample", &PythonAPI::resample,
             "Resample the raw data into a smoothed warped grid",
             py::arg("raw_data"), py::arg("num_mz") = 10,
             py::arg("num_rt") = 10, py::arg("smoothing_coef_mz") = 0.5,
             py::arg("smoothing_coef_rt") = 0.5)
        .def("resample_2", &PythonAPI::resample_2,
             "Resample the raw data into a smoothed warped grid",
             py::arg("raw_data"), py::arg("num_mz") = 10,
             py::arg("num_rt") = 10, py::arg("smoothing_coef_mz") = 0.5,
             py::arg("smoothing_coef_rt") = 0.5)
        .def("resample_3", &PythonAPI::resample_3,
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
