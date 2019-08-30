#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <thread>
#include <tuple>

#include "MIDAs/MIDAs.h"
#include "centroid/centroid.hpp"
#include "centroid/centroid_files.hpp"
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
#include "warp2d/warp2d.hpp"
#include "warp2d/warp2d_runners.hpp"

namespace py = pybind11;

namespace PythonAPI {
RawData::RawData read_mzxml(std::string file_name, double min_mz, double max_mz,
                            double min_rt, double max_rt,
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
            int64_t index = i + j * mesh.n;

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

    double theoretical_sigma_mz = RawData::fwhm_to_sigma(
        RawData::theoretical_fwhm(raw_data, mesh.bins_mz[local_max.i]));
    double theoretical_sigma_rt = RawData::fwhm_to_sigma(raw_data.fwhm_rt);
    {
        // Calculate the ROI for a given local max.
        double mz = peak.local_max_mz;
        double rt = peak.local_max_rt;

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
                if ((mz > peak.local_max_mz - theoretical_sigma_mz) &&
                    (mz < peak.local_max_mz + theoretical_sigma_mz) &&
                    (rt > peak.local_max_rt - theoretical_sigma_rt) &&
                    (rt < peak.local_max_rt + theoretical_sigma_rt)) {
                    peak.raw_roi_num_points_within_sigma++;
                }

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
                                       const Grid::Mesh &mesh,
                                       size_t max_peaks) {
    // Finding local maxima.
    auto local_max = find_local_max_idx(mesh);

    // The number of groups/threads is set to the maximum possible concurrency.
    uint64_t max_threads = std::thread::hardware_concurrency();

    // Split the points into different groups for concurrency.
    std::vector<std::vector<size_t>> groups =
        std::vector<std::vector<size_t>>(max_threads);
    for (size_t i = 0; i < local_max.size(); ++i) {
        size_t k = i % max_threads;
        groups[k].push_back(i);
    }

    std::vector<std::thread> threads(max_threads);
    std::vector<std::vector<Centroid::Peak>> peaks_array(max_threads);
    for (size_t i = 0; i < groups.size(); ++i) {
        threads[i] = std::thread(
            [&groups, &local_max, &peaks_array, &raw_data, &mesh, i]() {
                for (const auto &k : groups[i]) {
                    auto peak = build_peak(raw_data, mesh, local_max[k]);
                    // FIXME: Number of raw points within the theoretical sigma
                    // should be set by the user, with a sensible default. Same
                    // with the minimum number of rt scans per peak.
                    if (peak.raw_roi_num_points_within_sigma < 5) {
                        continue;
                    }
                    peaks_array[i].push_back(peak);
                }
            });
    }

    // Wait for the threads to finish.
    for (auto &thread : threads) {
        thread.join();
    }

    // Join peak groups.
    std::vector<Centroid::Peak> peaks;
    for (size_t i = 0; i < peaks_array.size(); ++i) {
        peaks.insert(end(peaks), begin(peaks_array[i]), end(peaks_array[i]));
    }

    // Sort the peaks by height.
    auto sort_peaks = [](const Centroid::Peak &p1,
                         const Centroid::Peak &p2) -> bool {
        return (p1.local_max_height > p2.local_max_height) ||
               (p1.local_max_height == p2.local_max_height);
    };
    std::stable_sort(peaks.begin(), peaks.end(), sort_peaks);

    // Update the peak ids.
    for (size_t i = 0; i < peaks.size(); ++i) {
        peaks[i].id = i;
    }

    // Return maximum amount of peaks.
    // TODO: Figure a way of performing max peaks when multiple threads are
    // in place without having to go through all of them. Perhaps an atomic
    // operation for increment counter?
    if (peaks.size() > max_peaks) {
        peaks.resize(max_peaks);
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
        return (p1.local_max_height > p2.local_max_height) ||
               (p1.local_max_height == p2.local_max_height);
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

void to_csv(const std::vector<Centroid::Peak> &peaks, std::string file_name) {
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

void write_raw_data(const RawData::RawData &raw_data, std::string file_name) {
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

RawData::RawData read_raw_data(std::string file_name) {
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

void write_mesh(const Grid::Mesh &mesh, std::string file_name) {
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

Grid::Mesh read_mesh(std::string file_name) {
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
                 std::string file_name) {
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

std::vector<Centroid::Peak> read_peaks(std::string file_name) {
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

void write_metamatch_clusters(
    const std::vector<MetaMatch::Cluster> &metamatch_clusters,
    std::string file_name) {
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

std::vector<MetaMatch::Cluster> read_metamatch_clusters(std::string file_name) {
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
                           std::string file_name) {
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

std::vector<MetaMatch::Peak> read_metamatch_peaks(std::string file_name) {
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

namespace IdentData {
struct SpectrumId {
    std::string id;
    bool pass_threshold;
    bool modifications;
    std::string sequence;
    std::string peptide_id;
    size_t charge_state;
    double theoretical_mz;
    double experimental_mz;
    double retention_time;
    int64_t rank;
};

struct CVParam {
    std::string name;
    std::string accession;
    std::string cv_ref;
    std::string value;
};

struct PeptideModification {
    double monoisotopic_mass_delta;
    double average_mass_delta;
    std::string residues;
    int64_t location;
    std::vector<CVParam> cv_params;
};

struct Peptide {
    std::string id;
    std::string sequence;
    std::vector<IdentData::PeptideModification> modifications;
};

struct DBSequence {
    std::string id;
    std::string value;
};

struct ProteinHypothesis {
    std::string db_sequence_id;
    bool pass_threshold;
    std::vector<std::string> spectrum_ids;
};

struct IdentData {
    std::vector<PythonAPI::IdentData::DBSequence> db_sequences;
    std::vector<PythonAPI::IdentData::Peptide> peptides;
    std::vector<PythonAPI::IdentData::SpectrumId> spectrum_ids;
    std::vector<PythonAPI::IdentData::ProteinHypothesis> protein_hypotheses;
};
}  // namespace IdentData

PythonAPI::IdentData::IdentData _read_mzidentml(std::istream &stream,
                                                bool threshold) {
    std::vector<PythonAPI::IdentData::SpectrumId> spectrum_ids;
    std::vector<PythonAPI::IdentData::Peptide> peptides;
    std::vector<PythonAPI::IdentData::DBSequence> db_sequences;
    std::vector<PythonAPI::IdentData::ProteinHypothesis> protein_hypotheses;

    // Find all Peptides and DBSequences in the file (SequenceCollection).
    while (stream.good()) {
        auto tag = XmlReader::read_tag(stream);
        if (!tag) {
            continue;
        }
        if (tag.value().name == "SequenceCollection" && tag.value().closed) {
            break;
        }
        if (tag.value().name == "Peptide") {
            auto peptide = PythonAPI::IdentData::Peptide{};
            auto attributes = tag.value().attributes;
            peptide.id = attributes["id"];
            while (stream.good()) {
                auto tag = XmlReader::read_tag(stream);
                if (tag.value().name == "Peptide" && tag.value().closed) {
                    peptides.push_back(peptide);
                    break;
                }
                if (tag.value().name == "PeptideSequence" &&
                    !tag.value().closed) {
                    auto data = XmlReader::read_data(stream);
                    if (!data) {
                        break;
                        // FIXME: Throw exception? Return nullopt?
                        // return std::nullopt;
                    }
                    peptide.sequence = data.value();
                }
                if (tag.value().name == "Modification" && !tag.value().closed) {
                    // Save modification info.
                    auto attributes = tag.value().attributes;
                    auto modification = IdentData::PeptideModification{};
                    if (attributes.find("monoisotopicMassDelta") !=
                        attributes.end()) {
                        modification.monoisotopic_mass_delta =
                            std::stod(attributes["monoisotopicMassDelta"]);
                    }
                    if (attributes.find("avgMassDelta") != attributes.end()) {
                        modification.average_mass_delta =
                            std::stod(attributes["avgMassDelta"]);
                    }
                    if (attributes.find("residues") != attributes.end()) {
                        modification.residues = attributes["residues"];
                    }
                    if (attributes.find("location") != attributes.end()) {
                        modification.location =
                            std::stoi(attributes["location"]);
                    }
                    // Find CVParams for this modification..
                    while (stream.good()) {
                        auto tag = XmlReader::read_tag(stream);
                        if (tag.value().name == "cvParam") {
                            auto cv_param = IdentData::CVParam{};
                            auto attributes = tag.value().attributes;
                            cv_param.name = attributes["name"];
                            cv_param.accession = attributes["accession"];
                            cv_param.cv_ref = attributes["cvRef"];
                            if (attributes.find("value") != attributes.end()) {
                                cv_param.value = attributes["value"];
                            }
                            modification.cv_params.push_back(cv_param);
                        }
                        if (tag.value().name == "Modification" &&
                            tag.value().closed) {
                            peptide.modifications.push_back(modification);
                            break;
                        }
                    }
                }
            }
        }
        if (tag.value().name == "DBSequence") {
            auto db_sequence = PythonAPI::IdentData::DBSequence{};
            auto attributes = tag.value().attributes;
            db_sequence.id = attributes["id"];
            while (stream.good()) {
                auto tag = XmlReader::read_tag(stream);
                if (tag.value().name == "DBSequence" && tag.value().closed) {
                    db_sequences.push_back(db_sequence);
                    break;
                }
                auto attributes = tag.value().attributes;
                if (tag.value().name == "cvParam" &&
                    attributes["accession"] == "MS:1001088") {
                    db_sequence.value = attributes["value"];
                }
            }
        }
    }
    // Find the PSMs for this data (SpectrumIdentificationResult).
    while (stream.good()) {
        auto tag = XmlReader::read_tag(stream);
        if (!tag) {
            continue;
        }
        if (tag.value().name == "SpectrumIdentificationList" &&
            tag.value().closed) {
            break;
        }
        if (tag.value().name != "SpectrumIdentificationResult") {
            continue;
        }
        auto spectrum_id = PythonAPI::IdentData::SpectrumId{};
        bool identification_item_found = false;
        while (stream.good()) {
            tag = XmlReader::read_tag(stream);
            auto attributes = tag.value().attributes;

            // Retention time.
            if (tag.value().name == "cvParam" &&
                attributes["accession"] == "MS:1000894") {
                spectrum_id.retention_time = std::stod(attributes["value"]);
            }

            // Identification item.
            if (tag.value().name == "SpectrumIdentificationItem" &&
                !tag.value().closed) {
                if (identification_item_found &&
                    std::stoi(attributes["rank"]) < spectrum_id.rank) {
                    continue;
                }
                spectrum_id.id = attributes["id"];
                spectrum_id.rank = std::stoi(attributes["rank"]);
                spectrum_id.pass_threshold =
                    attributes["passThreshold"] == "true";
                spectrum_id.peptide_id = attributes["peptide_ref"];
                spectrum_id.charge_state = std::stoi(attributes["chargeState"]);
                spectrum_id.theoretical_mz =
                    std::stod(attributes["calculatedMassToCharge"]);
                spectrum_id.experimental_mz =
                    std::stod(attributes["experimentalMassToCharge"]);
                identification_item_found = true;
            }

            if (tag.value().name == "SpectrumIdentificationResult" &&
                tag.value().closed) {
                if (!threshold || spectrum_id.pass_threshold) {
                    spectrum_ids.push_back(spectrum_id);
                }
                break;
            }
        }
    }
    // Find the protein groups for this data (ProteinDetectionList).
    while (stream.good()) {
        auto tag = XmlReader::read_tag(stream);
        if (!tag) {
            continue;
        }
        if (tag.value().name == "ProteinDetectionList" && tag.value().closed) {
            break;
        }
        if (tag.value().name == "ProteinDetectionHypothesis" &&
            !tag.value().closed) {
            auto protein_hypothesis = PythonAPI::IdentData::ProteinHypothesis{};
            auto attributes = tag.value().attributes;
            protein_hypothesis.db_sequence_id = attributes["dBSequence_ref"];
            protein_hypothesis.pass_threshold =
                attributes["passThreshold"] == "true";
            while (stream.good()) {
                tag = XmlReader::read_tag(stream);
                if (tag.value().name == "ProteinDetectionHypothesis" &&
                    tag.value().closed) {
                    if (!threshold || protein_hypothesis.pass_threshold) {
                        protein_hypotheses.push_back(protein_hypothesis);
                    }
                    break;
                }
                if (tag.value().name == "SpectrumIdentificationItemRef") {
                    auto attributes = tag.value().attributes;
                    protein_hypothesis.spectrum_ids.push_back(
                        attributes["spectrumIdentificationItem_ref"]);
                }
            }
        }
    }
    // Cross link peptide_id per SpectrumId to obtain the original sequence.
    for (auto &ident : spectrum_ids) {
        for (const auto &peptide : peptides) {
            if (peptide.id == ident.peptide_id) {
                ident.sequence = peptide.sequence;
                if (!peptide.modifications.empty()) {
                    ident.modifications = true;
                }
                break;
            }
        }
    }
    return {db_sequences, peptides, spectrum_ids, protein_hypotheses};
}  // namespace PythonAPI

PythonAPI::IdentData::IdentData read_mzidentml(std::string file_name,
                                               bool threshold) {
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
    return _read_mzidentml(stream, threshold);
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
    const std::vector<PythonAPI::IdentData::SpectrumId> &identifications,
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
    std::vector<double> sorted_rts;
    for (const auto &peak : sorted_peaks) {
        sorted_rts.push_back(peak.local_max_rt);
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

        // In case more than one peak is linked to the reference mz isotope, the
        // sequence with the less matching error should be selected. In order to
        // do so, the heights for each candidate must be normalized by the
        // reference isotope height.
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

                // Find the best matching candidate for the selected reference.
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
            // The weighted error for each normalized_candidate path can now be
            // evaluated.
            // NOTE: A potential alternative to an error function would be the
            // Euclidean distance between the expected theoretical peak and
            // candidates for that node.
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
        LinkedPeptide linked_peptide;
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

    py::class_<PythonAPI::IdentData::SpectrumId>(m, "SpectrumId")
        .def_readonly("id", &PythonAPI::IdentData::SpectrumId::id)
        .def_readonly("pass_threshold",
                      &PythonAPI::IdentData::SpectrumId::pass_threshold)
        .def_readonly("modifications",
                      &PythonAPI::IdentData::SpectrumId::modifications)
        .def_readonly("sequence", &PythonAPI::IdentData::SpectrumId::sequence)
        .def_readonly("peptide_id",
                      &PythonAPI::IdentData::SpectrumId::peptide_id)
        .def_readonly("charge_state",
                      &PythonAPI::IdentData::SpectrumId::charge_state)
        .def_readonly("theoretical_mz",
                      &PythonAPI::IdentData::SpectrumId::theoretical_mz)
        .def_readonly("experimental_mz",
                      &PythonAPI::IdentData::SpectrumId::experimental_mz)
        .def_readonly("retention_time",
                      &PythonAPI::IdentData::SpectrumId::retention_time)
        .def_readonly("rank", &PythonAPI::IdentData::SpectrumId::rank)
        .def("__repr__", [](const PythonAPI::IdentData::SpectrumId &s) {
            return s.sequence + "_" + std::to_string(s.charge_state);
        });

    py::class_<PythonAPI::IdentData::DBSequence>(m, "DBSequence")
        .def_readonly("id", &PythonAPI::IdentData::DBSequence::id)
        .def_readonly("value", &PythonAPI::IdentData::DBSequence::value);

    py::class_<PythonAPI::IdentData::PeptideModification>(m,
                                                          "PeptideModification")
        .def_readonly(
            "monoisotopic_mass_delta",
            &PythonAPI::IdentData::PeptideModification::monoisotopic_mass_delta)
        .def_readonly(
            "average_mass_delta",
            &PythonAPI::IdentData::PeptideModification::average_mass_delta)
        .def_readonly("residues",
                      &PythonAPI::IdentData::PeptideModification::residues)
        .def_readonly("location",
                      &PythonAPI::IdentData::PeptideModification::location);

    py::class_<PythonAPI::IdentData::Peptide>(m, "Peptide")
        .def_readonly("id", &PythonAPI::IdentData::Peptide::id)
        .def_readonly("sequence", &PythonAPI::IdentData::Peptide::sequence)
        .def_readonly("modifications",
                      &PythonAPI::IdentData::Peptide::modifications);

    py::class_<PythonAPI::IdentData::ProteinHypothesis>(m, "ProteinHypothesis")
        .def_readonly("db_sequence_id",
                      &PythonAPI::IdentData::ProteinHypothesis::db_sequence_id)
        .def_readonly("pass_threshold",
                      &PythonAPI::IdentData::ProteinHypothesis::pass_threshold)
        .def_readonly("spectrum_ids",
                      &PythonAPI::IdentData::ProteinHypothesis::spectrum_ids);

    py::class_<PythonAPI::IdentData::IdentData>(m, "IdentData")
        .def_readonly("db_sequences",
                      &PythonAPI::IdentData::IdentData::db_sequences)
        .def_readonly("peptides", &PythonAPI::IdentData::IdentData::peptides)
        .def_readonly("spectrum_ids",
                      &PythonAPI::IdentData::IdentData::spectrum_ids)
        .def_readonly("protein_hypotheses",
                      &PythonAPI::IdentData::IdentData::protein_hypotheses);

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
        .def("find_raw_points", &PythonAPI::find_raw_points,
             "Save the fitted peaks as a bpks file", py::arg("raw_data"),
             py::arg("min_mz"), py::arg("max_mz"), py::arg("min_rt"),
             py::arg("max_rt"))
        .def("find_peaks", &PythonAPI::find_peaks,
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
        .def("theoretical_isotopes_peptide",
             &PythonAPI::theoretical_isotopes_peptide,
             "Calculate the theoretical isotopic distribution of a peptide",
             py::arg("sequence"), py::arg("charge_state"),
             py::arg("min_perc") = 0.01)
        .def("read_mzidentml", &PythonAPI::read_mzidentml,
             "Read identification data from the given mzIdentML file ",
             py::arg("file_name"), py::arg("threshold") = true)
        .def("perform_metamatch", &PythonAPI::perform_metamatch,
             "Perform metamatch for peak matching", py::arg("input"),
             py::arg("radius_mz"), py::arg("radius_rt"), py::arg("fraction"))
        .def("link_identified_peptides", &PythonAPI::link_identified_peptides,
             "DEBUG", py::arg("peaks"), py::arg("identifications"),
             py::arg("tolerance_rt"), py::arg("minimum_isotope_perc"));
}
