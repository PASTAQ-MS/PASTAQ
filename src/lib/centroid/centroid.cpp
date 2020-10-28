#include <algorithm>
#include <thread>

#include "Eigen/Dense"

#include "centroid/centroid.hpp"
#include "utils/search.hpp"

#define PI 3.141592653589793238

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

    // Calculate the ROI for a given local max.
    double theoretical_sigma_mz = RawData::fwhm_to_sigma(
        RawData::theoretical_fwhm(raw_data, local_max.mz));
    double theoretical_sigma_rt = RawData::fwhm_to_sigma(raw_data.fwhm_rt);
    peak.roi_min_mz = peak.local_max_mz - 2 * theoretical_sigma_mz;
    peak.roi_max_mz = peak.local_max_mz + 2 * theoretical_sigma_mz;
    peak.roi_min_rt = peak.local_max_rt - 2 * theoretical_sigma_rt;
    peak.roi_max_rt = peak.local_max_rt + 2 * theoretical_sigma_rt;

    // Extract the raw data points for the ROI.
    auto raw_points =
        RawData::raw_points(raw_data, peak.roi_min_mz, peak.roi_max_mz,
                            peak.roi_min_rt, peak.roi_max_rt);
    if (raw_points.num_points == 0 || raw_points.num_scans < 3) {
        return std::nullopt;
    }

    {
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

    {
        // Solve the linearized 2D gaussian fitting problem `A * beta = c` with
        // weighted residuals.
        Eigen::MatrixXd A(5, 5);
        Eigen::VectorXd c(5);
        A = Eigen::MatrixXd::Zero(5, 5);
        c = Eigen::VectorXd::Zero(5);
        for (size_t i = 0; i < raw_points.num_points; ++i) {
            double mz = raw_points.mz[i] - local_max.mz;
            double rt = raw_points.rt[i] - local_max.rt;
            double intensity = raw_points.intensity[i];
            if (intensity <= 0) {
                continue;
            }

            double a = mz / theoretical_sigma_mz;
            double b = rt / theoretical_sigma_rt;
            double weight = intensity * std::exp(-0.5 * (a * a + b * b));
            double w_2 = weight * weight;

            A(0, 0) += w_2;
            A(0, 1) += w_2 * mz;
            A(0, 2) += w_2 * mz * mz;
            A(0, 3) += w_2 * rt;
            A(0, 4) += w_2 * rt * rt;

            A(1, 0) += w_2 * mz;
            A(1, 1) += w_2 * mz * mz;
            A(1, 2) += w_2 * mz * mz * mz;
            A(1, 3) += w_2 * rt * mz;
            A(1, 4) += w_2 * rt * rt * mz;

            A(2, 0) += w_2 * mz * mz;
            A(2, 1) += w_2 * mz * mz * mz;
            A(2, 2) += w_2 * mz * mz * mz * mz;
            A(2, 3) += w_2 * rt * mz * mz;
            A(2, 4) += w_2 * rt * rt * mz * mz;

            A(3, 0) += w_2 * rt;
            A(3, 1) += w_2 * mz * rt;
            A(3, 2) += w_2 * mz * mz * rt;
            A(3, 3) += w_2 * rt * rt;
            A(3, 4) += w_2 * rt * rt * rt;

            A(4, 0) += w_2 * rt * rt;
            A(4, 1) += w_2 * mz * rt * rt;
            A(4, 2) += w_2 * mz * mz * rt * rt;
            A(4, 3) += w_2 * rt * rt * rt;
            A(4, 4) += w_2 * rt * rt * rt * rt;

            c(0) += w_2 * std::log(intensity);
            c(1) += w_2 * std::log(intensity) * mz;
            c(2) += w_2 * std::log(intensity) * mz * mz;
            c(3) += w_2 * std::log(intensity) * rt;
            c(4) += w_2 * std::log(intensity) * rt * rt;
        }
        Eigen::VectorXd beta(5);
        beta = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(c);
        {
            double a = beta(0);
            double b = beta(1);
            double c = beta(2);
            double d = beta(3);
            double e = beta(4);
            if (std::isnan(a) || std::isnan(b) || std::isnan(c) ||
                std::isnan(d) || std::isnan(e) || std::isinf(a) ||
                std::isinf(b) || std::isinf(c) || std::isinf(d) ||
                std::isinf(e) || c >= 0 || e >= 0) {
                return std::nullopt;
            }
            double sigma_mz = std::sqrt(1 / (-2 * c));
            double mz = b / (-2 * c) + local_max.mz;
            double sigma_rt = std::sqrt(1 / (-2 * e));
            double rt = d / (-2 * e) + local_max.rt;
            double height =
                std::exp(a - ((b * b) / (4 * c)) - ((d * d) / (4 * e)));

            if (std::isnan(height) || std::isnan(mz) || std::isnan(sigma_mz) ||
                std::isnan(rt) || std::isnan(sigma_rt) || std::isinf(height) ||
                std::isinf(mz) || std::isinf(sigma_mz) || std::isinf(rt) ||
                std::isinf(sigma_rt) || height <= 0 || sigma_mz <= 0 ||
                sigma_rt <= 0 || mz <= 0 || rt <= 0) {
                return std::nullopt;
            }

            peak.fitted_height = height;
            peak.fitted_mz = mz;
            peak.fitted_rt = rt;
            peak.fitted_sigma_mz = sigma_mz;
            peak.fitted_sigma_rt = sigma_rt;
            peak.fitted_volume = height * sigma_mz * sigma_rt * 2.0 * PI;
        }
    }

    // Ensure peak quality.
    if (peak.raw_roi_sigma_mz <= 0 || peak.raw_roi_sigma_rt <= 0 ||
        peak.fitted_height > 2 * peak.raw_roi_max_height ||
        peak.fitted_mz < local_max.mz - 3 * theoretical_sigma_mz ||
        peak.fitted_mz > local_max.mz + 3 * theoretical_sigma_mz ||
        peak.fitted_rt < local_max.rt - 3 * theoretical_sigma_rt ||
        peak.fitted_rt > local_max.rt + 3 * theoretical_sigma_rt ||
        peak.fitted_sigma_mz <= theoretical_sigma_mz / 3 ||
        peak.fitted_sigma_rt <= theoretical_sigma_rt / 3 ||
        peak.fitted_sigma_mz >= theoretical_sigma_mz * 3 ||
        peak.fitted_sigma_rt >= theoretical_sigma_rt * 3) {
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
        return (p2.value < p1.value);
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
        return (p2.fitted_height < p1.fitted_height);
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
    double peak_a_mz = peak_a.fitted_mz;
    double peak_b_mz = peak_b.fitted_mz;
    double peak_a_rt = peak_a.fitted_rt + peak_a.rt_delta;
    double peak_b_rt = peak_b.fitted_rt + peak_b.rt_delta;
    // Early return if the peaks do not intersect in the +/-3 * sigma_mz/rt
    {
        double min_rt_a = peak_a_rt - 3 * peak_a.fitted_sigma_rt;
        double max_rt_a = peak_a_rt + 3 * peak_a.fitted_sigma_rt;
        double min_mz_a = peak_a_mz - 3 * peak_a.fitted_sigma_mz;
        double max_mz_a = peak_a_mz + 3 * peak_a.fitted_sigma_mz;
        double min_rt_b = peak_b_rt - 3 * peak_b.fitted_sigma_rt;
        double max_rt_b = peak_b_rt + 3 * peak_b.fitted_sigma_rt;
        double min_mz_b = peak_b_mz - 3 * peak_b.fitted_sigma_mz;
        double max_mz_b = peak_b_mz + 3 * peak_b.fitted_sigma_mz;

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

        return var_a * var_b * std::exp(0.5 * (a - b)) /
               std::sqrt(var_a + var_b);
    };

    auto rt_contrib = gaussian_contribution(
        peak_a_rt, peak_b_rt, peak_a.fitted_sigma_rt, peak_b.fitted_sigma_rt);
    auto mz_contrib = gaussian_contribution(
        peak_a_mz, peak_b_mz, peak_a.fitted_sigma_mz, peak_b.fitted_sigma_mz);

    return rt_contrib * mz_contrib * peak_a.fitted_height *
           peak_b.fitted_height;
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
