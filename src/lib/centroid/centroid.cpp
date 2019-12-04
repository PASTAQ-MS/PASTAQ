#include <algorithm>
#include <cmath>
#include <thread>

#include "centroid/centroid.hpp"
#include "utils/search.hpp"

std::vector<Centroid::LocalMax> Centroid::find_local_maxima(
    const Grid::Grid &grid) {
    std::vector<Centroid::LocalMax> points;
    // FIXME: This is performed in O(n^2), but using the divide and conquer
    // strategy we might achieve O(n * log(n)) or lower.
    // FIXME: Also, we should consider the corner case where neighbours are
    // exactly equal, both should be considered a local maxima and the average
    // of mz and rt should be reported.
    for (size_t j = 1; j < grid.m - 1; ++j) {
        for (size_t i = 1; i < grid.n - 1; ++i) {
            int64_t index = i + j * grid.n;

            // The definition of a local maxima in a 2D space might
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
            double value = grid.data[index];
            double right_value = grid.data[index + 1];
            double left_value = grid.data[index - 1];
            double top_value = grid.data[index - grid.n];
            double bottom_value = grid.data[index + grid.n];

            if ((value != 0) && (value > left_value) && (value > right_value) &&
                (value > top_value) && (value > bottom_value)) {
                points.push_back({grid.bins_mz[i], grid.bins_rt[j], value});
            }
        }
    }

    return points;
}

std::optional<Centroid::Peak> Centroid::build_peak(
    const RawData::RawData &raw_data, const LocalMax &local_max) {
    Centroid::Peak peak = {};
    peak.id = 0;
    peak.local_max_mz = local_max.mz;
    peak.local_max_rt = local_max.rt;
    peak.local_max_height = local_max.value;

    double theoretical_sigma_mz = RawData::fwhm_to_sigma(
        RawData::theoretical_fwhm(raw_data, local_max.mz));
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
        auto raw_points =
            RawData::raw_points(raw_data, peak.roi_min_mz, peak.roi_max_mz,
                                peak.roi_min_rt, peak.roi_max_rt);
        if (raw_points.num_points == 0) {
            return std::nullopt;
        }

        // Calculate the first 4 central moments for both mz/rt on the raw
        // data points using a 2 pass algorithm.
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
        for (size_t i = 0; i < raw_points.num_points; ++i) {
            double mz = raw_points.mz[i];
            double rt = raw_points.rt[i];
            double value = raw_points.intensity[i];
            if (value > max_value) {
                max_value = value;
            }
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
        if (weight_sum == 0) {
            return std::nullopt;
        }
        mz_mean /= weight_sum;
        rt_mean /= weight_sum;
        // Second pass.
        for (size_t i = 0; i < raw_points.num_points; ++i) {
            double mz = raw_points.mz[i];
            double rt = raw_points.rt[i];
            double value = raw_points.intensity[i];

            double mz_delta = mz - mz_mean;
            mz_m2 += value * std::pow(mz_delta, 2);
            mz_m3 += value * std::pow(mz_delta, 3);
            mz_m4 += value * std::pow(mz_delta, 4);

            double rt_delta = rt - rt_mean;
            rt_m2 += value * std::pow(rt_delta, 2);
            rt_m3 += value * std::pow(rt_delta, 3);
            rt_m4 += value * std::pow(rt_delta, 4);
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
        peak.raw_roi_num_points = raw_points.num_points;
        peak.raw_roi_num_scans = raw_points.num_scans;
        peak.raw_roi_total_intensity = weight_sum;
    }
    // FIXME: Number of raw points within the theoretical sigma
    // should be set by the user, with a sensible default. Same
    // with the minimum number of rt scans per peak.
    // Ensure peak quality.
    if (peak.raw_roi_num_points_within_sigma < 5 ||
        peak.raw_roi_num_scans < 3 || peak.raw_roi_sigma_mz == 0 ||
        peak.raw_roi_sigma_rt == 0) {
        return std::nullopt;
    }
    return peak;
}

std::vector<Centroid::Peak> Centroid::find_peaks_serial(
    const RawData::RawData &raw_data, const Grid::Grid &grid,
    size_t max_peaks) {
    // Finding local maxima.
    auto local_max = Centroid::find_local_maxima(grid);

    // Sort the local_maxima by value.
    auto sort_local_max = [](const Centroid::LocalMax &p1,
                             const Centroid::LocalMax &p2) -> bool {
        return (p1.value >= p2.value);
    };
    std::sort(local_max.begin(), local_max.end(), sort_local_max);

    std::vector<Centroid::Peak> peaks;
    for (const auto &max : local_max) {
        if (peaks.size() == max_peaks) {
            break;
        }
        auto peak = build_peak(raw_data, max);
        if (peak) {
            peaks.push_back(peak.value());
        }
    }

    // Update the peak ids.
    for (size_t i = 0; i < peaks.size(); ++i) {
        peaks[i].id = i;
    }

    return peaks;
}

std::vector<Centroid::Peak> Centroid::find_peaks_parallel(
    const RawData::RawData &raw_data, const Grid::Grid &grid, size_t max_peaks,
    size_t max_threads) {
    // Finding local maxima.
    auto local_max = Centroid::find_local_maxima(grid);

    // The number of groups/threads is set to the maximum possible concurrency.
    uint64_t num_threads = std::thread::hardware_concurrency();
    if (num_threads > max_threads) {
        num_threads = max_threads;
    }

    // Split the points into different groups for concurrency.
    std::vector<std::vector<size_t>> groups =
        std::vector<std::vector<size_t>>(num_threads);
    for (size_t i = 0; i < local_max.size(); ++i) {
        size_t k = i % num_threads;
        groups[k].push_back(i);
    }

    std::vector<std::thread> threads(num_threads);
    std::vector<std::vector<Centroid::Peak>> peaks_array(num_threads);
    for (size_t i = 0; i < groups.size(); ++i) {
        threads[i] = std::thread(
            [&groups, &local_max, &peaks_array, &raw_data, &grid, i]() {
                for (const auto &k : groups[i]) {
                    auto peak = build_peak(raw_data, local_max[k]);
                    if (peak) {
                        peaks_array[i].push_back(peak.value());
                    }
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
        return (p1.local_max_height >= p2.local_max_height);
    };
    std::sort(peaks.begin(), peaks.end(), sort_peaks);

    // Update the peak ids.
    for (size_t i = 0; i < peaks.size(); ++i) {
        peaks[i].id = i;
    }

    // Return maximum amount of peaks.
    if (peaks.size() > max_peaks) {
        peaks.resize(max_peaks);
    }

    return peaks;
}

double Centroid::peak_overlap(const Centroid::Peak &peak_a,
                              const Centroid::Peak &peak_b) {
    double peak_a_mz = peak_a.local_max_mz;
    double peak_b_mz = peak_b.local_max_mz;
    double peak_a_rt = peak_a.local_max_rt + peak_a.rt_delta;
    double peak_b_rt = peak_b.local_max_rt + peak_b.rt_delta;
    // Early return if the peaks do not intersect in the +/-3 * sigma_mz/rt
    {
        double min_rt_a = peak_a_rt - 3 * peak_a.raw_roi_sigma_rt;
        double max_rt_a = peak_a_rt + 3 * peak_a.raw_roi_sigma_rt;
        double min_mz_a = peak_a_mz - 3 * peak_a.raw_roi_sigma_mz;
        double max_mz_a = peak_a_mz + 3 * peak_a.raw_roi_sigma_mz;
        double min_rt_b = peak_b_rt - 3 * peak_b.raw_roi_sigma_rt;
        double max_rt_b = peak_b_rt + 3 * peak_b.raw_roi_sigma_rt;
        double min_mz_b = peak_b_mz - 3 * peak_b.raw_roi_sigma_mz;
        double max_mz_b = peak_b_mz + 3 * peak_b.raw_roi_sigma_mz;

        if (max_rt_a < min_rt_b || max_rt_b < min_rt_a || max_mz_a < min_mz_b ||
            max_mz_b < min_mz_a) {
            return 0;
        }
    }

    // Calculate the gaussian contribution of the overlap between two points in
    // one dimension.
    auto gaussian_contribution = [](double x_a, double x_b, double sigma_a,
                                    double sigma_b) -> double {
        double var_a = std::pow(sigma_a, 2);
        double var_b = std::pow(sigma_b, 2);

        double a = (var_a + var_b) / (var_a * var_b) *
                   std::pow((x_a * var_b + x_b * var_a) / (var_a + var_b), 2);
        double b = (x_a * x_a) / var_a + (x_b * x_b) / var_b;

        return std::exp(0.5 * (a - b)) / std::sqrt(var_a + var_b);
    };

    auto rt_contrib = gaussian_contribution(
        peak_a_rt, peak_b_rt, peak_a.raw_roi_sigma_rt, peak_b.raw_roi_sigma_rt);
    auto mz_contrib = gaussian_contribution(
        peak_a_mz, peak_b_mz, peak_a.raw_roi_sigma_mz, peak_b.raw_roi_sigma_mz);

    return rt_contrib * mz_contrib * peak_a.local_max_height *
           peak_b.local_max_height;
}

double Centroid::cumulative_overlap(const std::vector<Centroid::Peak> &set_a,
                                    const std::vector<Centroid::Peak> &set_b) {
    double total_overlap = 0;
    for (const auto &peak_a : set_a) {
        for (const auto &peak_b : set_b) {
            total_overlap += Centroid::peak_overlap(peak_a, peak_b);
        }
    }
    return total_overlap;
}
