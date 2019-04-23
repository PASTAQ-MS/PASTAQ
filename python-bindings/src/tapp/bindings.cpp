#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <tuple>

#include "centroid/centroid.hpp"
#include "centroid/centroid_files.hpp"
#include "grid/grid.hpp"
#include "grid/grid_files.hpp"
#include "grid/raw_data.hpp"
#include "grid/xml_reader.hpp"
#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "utils/search.hpp"
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

    if (min_rt < raw_data.min_rt) {
        min_rt = raw_data.min_rt;
    }
    if (max_rt > raw_data.max_rt) {
        max_rt = raw_data.max_rt;
    }
    size_t min_j = Search::lower_bound(raw_data.retention_times, min_rt);
    size_t max_j = scans.size();
    if (scans[min_j].retention_time < min_rt) {
        ++min_j;
    }

    for (size_t j = min_j; j < max_j; ++j) {
        const auto &scan = scans[j];
        if (scan.num_points == 0) {
            continue;
        }
        if (scan.retention_time > max_rt) {
            break;
        }

        size_t min_i = Search::lower_bound(scan.mz, min_mz);
        size_t max_i = scan.num_points;
        if (scan.mz[min_i] < min_mz) {
            ++min_i;
        }
        for (size_t i = min_i; i < max_i; ++i) {
            if (scan.mz[i] > max_mz) {
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

struct MeshIndex {
    size_t i;
    size_t j;
};

std::vector<MeshIndex> find_local_max_idx(const Grid::Mesh &mesh) {
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
            double value = mesh.data[index];
            double right_value = mesh.data[index + 1];
            double left_value = mesh.data[index - 1];
            double top_value = mesh.data[index - mesh.n];
            double bottom_value = mesh.data[index + mesh.n];

            if ((value != 0) && (value > left_value) && (value > right_value) &&
                (value > top_value) && (value > bottom_value)) {
                points.push_back({i, j});
            }
        }
    }

    auto sort_by_value = [&mesh](const MeshIndex &p1,
                                 const MeshIndex &p2) -> bool {
        return (mesh.data[p1.i + p1.j * mesh.n] >
                mesh.data[p2.i + p2.j * mesh.n]);
    };
    std::stable_sort(points.begin(), points.end(), sort_by_value);

    return points;
}

void explore_peak_slope(uint64_t i, uint64_t j, double previous_value,
                        const Grid::Mesh &mesh,
                        std::vector<MeshIndex> &points) {
    // Check that the point has not being already included.
    for (const auto &point : points) {
        if (point.i == i && point.j == j) {
            return;
        }
    }

    double value = mesh.data[i + j * mesh.n];
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

Centroid::Peak build_peak(const RawData::RawData &raw_data,
                          const Grid::Mesh &mesh, const MeshIndex &local_max) {
    Centroid::Peak peak = {};
    peak.id = 0;
    peak.local_max_mz = mesh.bins_mz[local_max.i];
    peak.local_max_rt = mesh.bins_rt[local_max.j];
    peak.local_max_height = mesh.data[local_max.i + local_max.j * mesh.n];

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
            boundary_sum += mesh.data[point.i + point.j * mesh.n];
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
            double value = mesh.data[point.i + point.j * mesh.n];

            height_sum += value;
            x_sum += value * mz;
            y_sum += value * rt;
            x_sig += value * mz * mz;
            y_sig += value * rt * rt;
        }
        peak.slope_descent_mean_mz = x_sum / height_sum;
        peak.slope_descent_mean_rt = y_sum / height_sum;
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

        double theoretical_sigma_mz = RawData::fwhm_to_sigma(
            RawData::theoretical_fwhm(raw_data, mesh.bins_mz[local_max.i]));
        double theoretical_sigma_rt = RawData::fwhm_to_sigma(raw_data.fwhm_rt);

        peak.roi_min_mz = mz - 3 * theoretical_sigma_mz;
        peak.roi_max_mz = mz + 3 * theoretical_sigma_mz;
        peak.roi_min_rt = rt - 3 * theoretical_sigma_rt;
        peak.roi_max_rt = rt + 3 * theoretical_sigma_rt;
    }

    {
        const auto &scans = raw_data.scans;
        // FIXME: Make nan instead?
        if (scans.size() == 0) {
            std::ostringstream error_stream;
            error_stream << "the given raw_data is empty";
            throw std::invalid_argument(error_stream.str());
        }

        size_t min_j =
            Search::lower_bound(raw_data.retention_times, peak.roi_min_rt);
        size_t max_j = scans.size();
        if (scans[min_j].retention_time < peak.roi_min_rt) {
            ++min_j;
        }

        // Calculate the first 4 central moments for both mz/rt on the raw
        // data points using a 2 pass algorithm.
        size_t num_scans = 0;
        size_t num_points = 0;
        double max_value = 0;
        double mz_mean = 0;
        double mz_m2 = 0;
        double mz_m3 = 0;
        double mz_m4 = 0;
        double rt_mean = 0;
        double rt_m2 = 0;
        double rt_m3 = 0;
        double rt_m4 = 0;
        double weight_sum = 0;
        // First pass.
        for (size_t j = min_j; j < max_j; ++j) {
            const auto &scan = scans[j];
            if (scan.retention_time > peak.roi_max_rt) {
                break;
            }
            if (scan.num_points == 0) {
                continue;
            }

            size_t min_i = Search::lower_bound(scan.mz, peak.roi_min_mz);
            size_t max_i = scan.num_points;
            if (scan.mz[min_i] < peak.roi_min_mz) {
                ++min_i;
            }
            bool scan_not_empty = false;
            for (size_t i = min_i; i < max_i; ++i) {
                if (scan.mz[i] > peak.roi_max_mz) {
                    break;
                }
                double mz = scan.mz[i];
                double rt = scan.retention_time;
                double value = scan.intensity[i];
                if (value > max_value) {
                    max_value = value;
                }
                scan_not_empty = true;
                ++num_points;

                weight_sum += value;
                mz_mean += value * mz;
                rt_mean += value * rt;
            }
            if (scan_not_empty) {
                ++num_scans;
            }
        }
        mz_mean /= weight_sum;
        rt_mean /= weight_sum;

        // Second pass.
        for (size_t j = min_j; j < max_j; ++j) {
            const auto &scan = scans[j];
            if (scan.retention_time > peak.roi_max_rt) {
                break;
            }
            if (scan.num_points == 0) {
                continue;
            }

            size_t min_i = Search::lower_bound(scan.mz, peak.roi_min_mz);
            size_t max_i = scan.num_points;
            if (scan.mz[min_i] < peak.roi_min_mz) {
                ++min_i;
            }
            for (size_t i = min_i; i < max_i; ++i) {
                if (scan.mz[i] > peak.roi_max_mz) {
                    break;
                }
                double mz = scan.mz[i];
                double rt = scan.retention_time;
                double value = scan.intensity[i];

                double mz_delta = mz - mz_mean;
                mz_m2 += value * std::pow(mz_delta, 2);
                mz_m3 += value * std::pow(mz_delta, 3);
                mz_m4 += value * std::pow(mz_delta, 4);

                double rt_delta = rt - rt_mean;
                rt_m2 += value * std::pow(rt_delta, 2);
                rt_m3 += value * std::pow(rt_delta, 3);
                rt_m4 += value * std::pow(rt_delta, 4);
            }
        }
        mz_m2 /= weight_sum;
        mz_m3 /= weight_sum;
        mz_m4 /= weight_sum;
        rt_m2 /= weight_sum;
        rt_m3 /= weight_sum;
        rt_m4 /= weight_sum;

        // Update the peak data structure.
        peak.raw_roi_mean_mz = mz_mean;
        peak.raw_roi_sigma_mz = std::sqrt(mz_m2);
        peak.raw_roi_skewness_mz = mz_m3 / std::pow(mz_m2, 1.5);
        peak.raw_roi_kurtosis_mz = mz_m4 / std::pow(mz_m2, 2);
        peak.raw_roi_mean_rt = rt_mean;
        peak.raw_roi_sigma_rt = std::sqrt(rt_m2);
        peak.raw_roi_skewness_rt = rt_m3 / std::pow(rt_m2, 1.5);
        peak.raw_roi_kurtosis_rt = rt_m4 / std::pow(rt_m2, 2);
        peak.raw_roi_max_height = max_value;
        peak.raw_roi_total_intensity = weight_sum;
        peak.raw_roi_num_points = num_points;
        peak.raw_roi_num_scans = num_scans;
        peak.raw_roi_total_intensity = weight_sum;

        // FIXME: Make nan instead?
        // if (raw_points.num_points == 0) {
        // std::ostringstream error_stream;
        // error_stream << "couldn't find raw_data points on the given ROI";
        // throw std::invalid_argument(error_stream.str());
        //}
    }

    return peak;
}

std::vector<Centroid::Peak> find_peaks(const RawData::RawData &raw_data,
                                       const Grid::Mesh &mesh) {
    auto local_max = find_local_max_idx(mesh);
    std::vector<Centroid::Peak> peaks;
    size_t i = 0;
    for (const auto &lm : local_max) {
        auto peak = build_peak(raw_data, mesh, lm);
        peak.id = i;
        peaks.push_back(peak);
        ++i;
    }
    return peaks;
}

std::vector<std::vector<Centroid::Peak>> warp_peaks(
    const std::vector<std::vector<Centroid::Peak>> &all_peaks,
    size_t reference_index, int64_t slack, int64_t window_size,
    int64_t num_points, double rt_expand_factor, int64_t peaks_per_window) {
    // TODO(alex): Validate the parameters and throw an error if appropriate.
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
            Warp2D::Runners::Serial::run(reference_peaks, peaks, parameters);
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
        return (p1.local_max_height > p2.local_max_height) ||
               (p1.local_max_height == p2.local_max_height);
    };
    std::stable_sort(peak_list_a.begin(), peak_list_a.end(), sort_peaks);
    std::stable_sort(peak_list_b.begin(), peak_list_b.end(), sort_peaks);
    peak_list_a.resize(n_peaks);
    peak_list_b.resize(n_peaks);
    SimilarityResults results = {};
    results.self_a = Warp2D::similarity_2D(peak_list_a, peak_list_a);
    results.self_b = Warp2D::similarity_2D(peak_list_b, peak_list_b);
    results.overlap = Warp2D::similarity_2D(peak_list_a, peak_list_b);
    // Overlap / (GeometricMean(self_a, self_b))
    results.geometric_ratio =
        results.overlap / std::sqrt(results.self_a * results.self_b);
    // Harmonic mean of the ratios between self_similarity/overlap_similarity
    results.mean_ratio =
        2 * results.overlap / (results.self_a + results.self_b);
    return results;
}

struct PeakList {
    std::vector<Centroid::Peak> peaks;
    std::string file_name;   // NOTE: Should this be on the raw_data instead?
    std::string class_name;  // NOTE: Should this be on the raw_data instead?
    std::shared_ptr<RawData::RawData> raw_data;
};

struct WarpingTimeMap {
    std::vector<double> rt_start;
    std::vector<double> rt_end;
    std::vector<double> warped_rt_start;
    std::vector<double> warped_rt_end;
    size_t n_warping_segments;
};

struct PeakLists {
    std::vector<PeakList> peak_lists;
    std::shared_ptr<WarpingTimeMap> warping_time_map;
};

struct MetaPeak {
    // TODO: ...
};

struct MetaPeaks {
    // TODO: ...
};

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

    py::class_<Grid::Mesh>(m, "Mesh")
        .def_readonly("n", &Grid::Mesh::n)
        .def_readonly("m", &Grid::Mesh::m)
        .def_readonly("data", &Grid::Mesh::data)
        .def_readonly("bins_mz", &Grid::Mesh::bins_mz)
        .def_readonly("bins_rt", &Grid::Mesh::bins_rt);

    py::class_<PythonAPI::RawPoints>(m, "RawPoints")
        .def_readonly("rt", &PythonAPI::RawPoints::rt)
        .def_readonly("mz", &PythonAPI::RawPoints::mz)
        .def_readonly("intensity", &PythonAPI::RawPoints::intensity);

    py::class_<Centroid::Peak>(m, "Peak")
        .def_readonly("id", &Centroid::Peak::id)
        .def_readonly("local_max_mz", &Centroid::Peak::local_max_mz)
        .def_readonly("local_max_rt", &Centroid::Peak::local_max_rt)
        .def_readonly("local_max_height", &Centroid::Peak::local_max_height)
        .def_readonly("slope_descent_mean_mz",
                      &Centroid::Peak::slope_descent_mean_mz)
        .def_readonly("slope_descent_mean_rt",
                      &Centroid::Peak::slope_descent_mean_rt)
        .def_readonly("slope_descent_sigma_mz",
                      &Centroid::Peak::slope_descent_sigma_mz)
        .def_readonly("slope_descent_sigma_rt",
                      &Centroid::Peak::slope_descent_sigma_rt)
        .def_readonly("slope_descent_total_intensity",
                      &Centroid::Peak::slope_descent_total_intensity)
        .def_readonly("slope_descent_border_background",
                      &Centroid::Peak::slope_descent_border_background)
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
             py::arg("method") = "sum");

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

    // Functions.
    m.def("read_mzxml", &PythonAPI::read_mzxml,
          "Read raw data from the given mzXML file ", py::arg("file_name"),
          py::arg("min_mz") = -1.0, py::arg("max_mz") = -1.0,
          py::arg("min_rt") = -1.0, py::arg("max_rt") = -1.0,
          py::arg("instrument_type") = "", py::arg("resolution_ms1"),
          py::arg("resolution_msn"), py::arg("reference_mz"),
          py::arg("fwhm_rt"), py::arg("polarity") = "")
        .def("theoretical_fwhm", &RawData::theoretical_fwhm,
             "Calculate the theoretical width of the peak at the given m/z for "
             "the given raw file",
             py::arg("raw_data"), py::arg("mz"))
        .def("resample", &Grid::resample,
             "Resample the raw data into a smoothed warped grid",
             py::arg("raw_data"), py::arg("num_mz") = 10,
             py::arg("num_rt") = 10, py::arg("smoothing_coef_mz") = 0.5,
             py::arg("smoothing_coef_rt") = 0.5)
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
             py::arg("rt_expand_factor"), py::arg("peaks_per_window"))
        .def("find_similarity", &PythonAPI::find_similarity,
             "Find the similarity between two peak lists",
             py::arg("peak_list_a"), py::arg("peak_list_b"),
             py::arg("n_peaks"));
}
