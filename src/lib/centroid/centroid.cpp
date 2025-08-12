#include <algorithm>
#include <thread>
#include <fstream>
#include <iostream>

#include "Eigen/Core"
#include "Eigen/Dense"

#include "centroid/centroid.hpp"
// #include "utils/search.hpp"

#define PI 3.141592653589793238

/* bool triggerPeakLog{true};

void log_itg(const std::string& msg) {
    const char* log_path = "/home/phorvatovich/development/pastaqTesting/PASTAQ/log.txt";
    // Open in append mode so previous logs are preserved
    if (msg == "reset") {
        // Delete the log file and start fresh
        std::remove(log_path);
    }

    std::ofstream log_file(log_path, std::ios::app);

    if (!log_file) {
        std::cerr << "Error: Could not open log file." << std::endl;
        return;
    }

    log_file << msg << std::endl;

    // Optional: flush explicitly
    log_file.flush();
} */

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
            if (raw_data.get_failed_peaks){
                peak.fit_failure_code |= 1 << 0;
            } else {
                return std::nullopt;
            }
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

        double sigma_mz{0};
        double mz{0};
        double sigma_rt{0};
        double rt{0};
        double height{0};
        /* std::string log_msg = std::string("Raw data type profile (0) or centroid (1): ") +
                std::to_string(raw_data.centroid) + "\n\n";
        log_itg(log_msg); */

        if (!(raw_data.centroid)) {
            /* std::string log_msg = std::string("Profile data peak picking.") + "\n\n";
            log_itg(log_msg); */

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
                    if (raw_data.get_failed_peaks){
                        peak.fit_failure_code |= static_cast<uint64_t>(std::isnan(a)) << 1;
                        peak.fit_failure_code  |= static_cast<uint64_t>(std::isnan(b)) << 2;
                        peak.fit_failure_code  |= static_cast<uint64_t>(std::isnan(c)) << 3;
                        peak.fit_failure_code  |= static_cast<uint64_t>(std::isnan(d)) << 4;
                        peak.fit_failure_code  |= static_cast<uint64_t>(std::isnan(e)) << 5;
                    
                        peak.fit_failure_code  |= static_cast<uint64_t>(std::isinf(a)) << 6;
                        peak.fit_failure_code  |= static_cast<uint64_t>(std::isinf(b)) << 7;
                        peak.fit_failure_code  |= static_cast<uint64_t>(std::isinf(c)) << 8;
                        peak.fit_failure_code  |= static_cast<uint64_t>(std::isinf(d)) << 9;
                        peak.fit_failure_code  |= static_cast<uint64_t>(std::isinf(e)) << 10;
                    
                        peak.fit_failure_code  |= static_cast<uint64_t>(c >= 0) << 11;
                        peak.fit_failure_code  |= static_cast<uint64_t>(e >= 0) << 12;
                    } else {
                        return std::nullopt;
                    }
                }
                sigma_mz = std::sqrt(1 / (-2 * c));
                mz = b / (-2 * c) + local_max.mz;            //it use 1/(-2c) in order to use the original value from the solved equation
                sigma_rt = std::sqrt(1 / (-2 * e));
                rt = d / (-2 * e) + local_max.rt;            // it use 1/(-2e) in order to use the original value from the solved equation
                height = std::exp(a - ((b * b) / (4 * c)) - ((d * d) / (4 * e)));

                /* log_msg = std::string("derived parameters. sigma mz: " + std::to_string(sigma_mz) + ", mz: " + std::to_string(mz) + ", sigma rt: " +
                    std::to_string(sigma_rt) + ", rt: " + std::to_string(rt) + ", height: " +
                    std::to_string(height)) + "\n" +
                    "a: " + std::to_string(a) + ", b: " + std::to_string(b) + ", c: " +
                    std::to_string(c) + ", d: " + std::to_string(d) + ", e: " + std::to_string(e) + "\n\n";

                if (local_max.mz>612.95&&local_max.mz<612.99&&local_max.rt>3635&&local_max.rt<3645) {
                    log_msg += "Roi range:\n";
                    log_msg += "mzRoiMin: " + std::to_string(peak.roi_min_mz) +
                        ", mzRoiMax: " + std::to_string(peak.roi_max_mz) +
                        ", rtRoiMin: " + std::to_string(peak.roi_min_rt) +
                        ", rtRoiMax: " + std::to_string(peak.roi_max_rt) + "\n\n";
                }
                
                if (local_max.mz>612.95&&local_max.mz<612.99&&local_max.rt>3635&&local_max.rt<3645) {
                    log_msg += "Raw data points:\n";
                    for (size_t i = 0; i < raw_points.num_points; ++i) {
                        log_msg += "mz: " + std::to_string(raw_points.mz[i]) +
                        " rt: " + std::to_string(raw_points.rt[i]) + " intensity: " +
                        std::to_string(raw_points.intensity[i]) + "\n";
                    }

                    log_msg += "\n";
                }
                log_itg(log_msg); */
            }
        }
        else {
            /* std::string log_msg = std::string("Centroid data peak picking.") + "\n\n";
            log_itg(log_msg); */
            
            // If the raw data is centroided, we can use simpler approach
            // to solve the linear system by fixing sigmamz and umz
            Eigen::VectorXd beta(3);
            Eigen::VectorXi selectedCols(3);
            selectedCols << 0, 3, 4; // index vector for selected columns
            Eigen::VectorXi fixedCols(2);
            fixedCols << 1, 2; // index vector for not-selected columns
            Eigen::MatrixXd B(5, 3);
            B = A(Eigen::indexing::all, selectedCols);
            Eigen::VectorXd e(2);
            e(0) = 0; //local_max.mz/(theoretical_sigma_mz*theoretical_sigma_mz); // these are the fixed constant the peak location in mz and sigma_mz
            e(1) = -1/(2*theoretical_sigma_mz*theoretical_sigma_mz);
            Eigen::VectorXd d = c - A(Eigen::indexing::all, fixedCols) * e;
            beta = B.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(d);

            /* // if (!triggerPeakLog) {
                log_msg = std::string("Beta values: ") + std::to_string(beta(0)) + ", " + std::to_string(beta(1)) + ", " + std::to_string(beta(2)) + "\n" +
                    "d vector: " + std::to_string(d(0)) + ", " + std::to_string(d(1)) + ", " + std::to_string(d(2)) + "\n" +
                    "e vector: " + std::to_string(e(0)) + ", " + std::to_string(e(1)) + "\n\n";
                log_itg(log_msg);
                // triggerPeakLog = false;
            // } */
            {
                double a = beta(0);
                double b = beta(1);
                double c = beta(2);
                if (std::isnan(a) || std::isnan(b) || std::isnan(c) ||
                    std::isinf(a) || std::isinf(b) || std::isinf(c) ||
                    c >= 0) {
                        if (raw_data.get_failed_peaks){
                            peak.fit_failure_code |= static_cast<uint64_t>(std::isnan(a)) << 1;
                            peak.fit_failure_code  |= static_cast<uint64_t>(std::isnan(b)) << 2;
                            peak.fit_failure_code  |= static_cast<uint64_t>(std::isnan(c)) << 3;
                        
                            peak.fit_failure_code  |= static_cast<uint64_t>(std::isinf(a)) << 6;
                            peak.fit_failure_code  |= static_cast<uint64_t>(std::isinf(b)) << 7;
                            peak.fit_failure_code  |= static_cast<uint64_t>(std::isinf(c)) << 8;
                        
                            peak.fit_failure_code  |= static_cast<uint64_t>(c >= 0) << 11;
                        } else {
                            return std::nullopt;
                        }
                }
                sigma_mz = theoretical_sigma_mz;
                mz = local_max.mz;
                sigma_rt = std::sqrt(1 / (-2 * c));
                rt = b / (-2 * c) + local_max.rt;
                height = std::exp(a - ((b * b) / (4 * c))); //+ (0.5 * (local_max.mz * local_max.mz) / (theoretical_sigma_mz * theoretical_sigma_mz)) 

                /* log_msg = std::string("derived parameters. sigma mz: " + std::to_string(sigma_mz) + ", mz: " + std::to_string(mz) + ", sigma rt: " +
                    std::to_string(sigma_rt) + ", rt: " + std::to_string(rt) + ", height: " +
                    std::to_string(height)) + "\n" +
                    "a: " + std::to_string(a) + ", b: " + std::to_string(b) + ", c: " +
                    std::to_string(c) + "\n" +
                    "height elements. First:  " + std::to_string(a) + ", Second: " + std::to_string((0.5 * (local_max.mz * local_max.mz) / (theoretical_sigma_mz * theoretical_sigma_mz))) + ", Third: " +
                    std::to_string(-((b * b) / (4 * c))) + ", Sum: " + std::to_string(a + (0.5 * (local_max.mz * local_max.mz) / (theoretical_sigma_mz * theoretical_sigma_mz)) - ((b * b) / (4 * c))) + "\n\n";
                
                if (local_max.mz>612.95&&local_max.mz<612.99&&local_max.rt>3635&&local_max.rt<3645) {
                    int num_points_per_lines = 10;
                    log_msg += "Raw data points:\n";
                    if (raw_points.num_points > num_points_per_lines) {
                        for (size_t i = 0; i < raw_points.num_points; i += num_points_per_lines) {
                            for (size_t j = 0; j < num_points_per_lines && i + j < raw_points.num_points; ++j) {
                                log_msg += "mz: " + std::to_string(raw_points.mz[i + j]) +
                                " rt: " + std::to_string(raw_points.rt[i + j]) + " intensity: " +
                                std::to_string(raw_points.intensity[i + j]) + "\n";
                            }
                            log_msg += "\n";
                        }
                    } else {
                        for (size_t i = 0; i < raw_points.num_points; ++i) {
                            log_msg += "mz: " + std::to_string(raw_points.mz[i]) +
                            " rt: " + std::to_string(raw_points.rt[i]) + " intensity: " +
                            std::to_string(raw_points.intensity[i]) + "\n";
                        }
                    }
                    log_msg += "\n";
                }

                log_itg(log_msg); */
            }
        }

        if (std::isnan(height) || std::isnan(mz) || std::isnan(sigma_mz) ||
        std::isnan(rt) || std::isnan(sigma_rt) || std::isinf(height) ||
        std::isinf(mz) || std::isinf(sigma_mz) || std::isinf(rt) ||
        std::isinf(sigma_rt) || height <= 0 || sigma_mz <= 0 ||
        sigma_rt <= 0 || mz <= 0 || rt <= 0) {
            if (raw_data.get_failed_peaks){
                peak.fit_failure_code  |= static_cast<uint64_t>(std::isnan(height)) << 13;
                peak.fit_failure_code |= static_cast<uint64_t>(std::isnan(mz)) << 14;
                peak.fit_failure_code |= static_cast<uint64_t>(std::isnan(sigma_mz)) << 15;
                peak.fit_failure_code |= static_cast<uint64_t>(std::isnan(rt)) << 16;
                peak.fit_failure_code |= static_cast<uint64_t>(std::isnan(sigma_rt)) << 17;
                peak.fit_failure_code |= static_cast<uint64_t>(std::isinf(height)) << 18;
                peak.fit_failure_code |= static_cast<uint64_t>(std::isinf(mz)) << 19;
                peak.fit_failure_code |= static_cast<uint64_t>(std::isinf(sigma_mz)) << 20;
                peak.fit_failure_code |= static_cast<uint64_t>(std::isinf(rt)) << 21;
                peak.fit_failure_code |= static_cast<uint64_t>(std::isinf(sigma_rt)) << 22;
                peak.fit_failure_code |= static_cast<uint64_t>(height <= 0) << 23;
                peak.fit_failure_code |= static_cast<uint64_t>(sigma_mz <= 0) << 24;
                peak.fit_failure_code |= static_cast<uint64_t>(sigma_rt <= 0) << 25;
                peak.fit_failure_code |= static_cast<uint64_t>(mz <= 0) << 26;
                peak.fit_failure_code |= static_cast<uint64_t>(rt <= 0) << 27;
                } else {
                    return std::nullopt;
            }
    }

        peak.fitted_height = height;
        peak.fitted_mz = mz;
        peak.fitted_rt = rt;
        peak.fitted_sigma_mz = sigma_mz;
        peak.fitted_sigma_rt = sigma_rt;
        peak.fitted_volume = height * sigma_mz * sigma_rt * 2.0 * PI;
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
        if (raw_data.get_failed_peaks){
            peak.fit_failure_code  |= static_cast<uint64_t>(peak.raw_roi_sigma_mz <= 0) << 28;
            peak.fit_failure_code  |= static_cast<uint64_t>(peak.raw_roi_sigma_rt <= 0) << 29;
            peak.fit_failure_code  |= static_cast<uint64_t>(peak.fitted_height > 2 * peak.raw_roi_max_height) << 30;
            peak.fit_failure_code  |= static_cast<uint64_t>(peak.fitted_mz < local_max.mz - 3 * theoretical_sigma_mz) << 31;
            peak.fit_failure_code  |= static_cast<uint64_t>(peak.fitted_mz > local_max.mz + 3 * theoretical_sigma_mz) << 32;
            peak.fit_failure_code  |= static_cast<uint64_t>(peak.fitted_rt < local_max.rt - 3 * theoretical_sigma_rt) << 33;
            peak.fit_failure_code  |= static_cast<uint64_t>(peak.fitted_rt > local_max.rt + 3 * theoretical_sigma_rt) << 34;
            peak.fit_failure_code  |= static_cast<uint64_t>(peak.fitted_sigma_mz <= theoretical_sigma_mz / 3) << 35;
            peak.fit_failure_code  |= static_cast<uint64_t>(peak.fitted_sigma_rt <= theoretical_sigma_rt / 3) << 36;
            peak.fit_failure_code  |= static_cast<uint64_t>(peak.fitted_sigma_mz >= theoretical_sigma_mz * 3) << 37;
            peak.fit_failure_code  |= static_cast<uint64_t>(peak.fitted_sigma_rt >= theoretical_sigma_rt * 3) << 38;
        } else {
            return std::nullopt;
        }
    }
    return peak;
}

const std::vector<std::string> Centroid::Peak::error_messages = {
    "Weight sum is zero",
    "NaN detected in parameter 'a'",
    "NaN detected in parameter 'b'",
    "NaN detected in parameter 'c'",
    "NaN detected in parameter 'd'",
    "NaN detected in parameter 'e'",
    "Infinity detected in parameter 'a'",
    "Infinity detected in parameter 'b'",
    "Infinity detected in parameter 'c'",
    "Infinity detected in parameter 'd'",
    "Infinity detected in parameter 'e'",
    "Parameter 'c' is non-negative",
    "Parameter 'e' is non-negative",
    "NaN detected in fitted height",
    "NaN detected in fitted mz",
    "NaN detected in fitted sigma_mz",
    "NaN detected in fitted rt",
    "NaN detected in fitted sigma_rt",
    "Infinity detected in fitted height",
    "Infinity detected in fitted mz",
    "Infinity detected in fitted sigma_mz",
    "Infinity detected in fitted rt",
    "Infinity detected in fitted sigma_rt",
    "Fitted height is non-positive",
    "Fitted sigma_mz is non-positive",
    "Fitted sigma_rt is non-positive",
    "Fitted mz is non-positive",
    "Fitted rt is non-positive",
    "Raw ROI sigma_mz is non-positive",
    "Raw ROI sigma_rt is non-positive",
    "Fitted height exceeds twice the raw ROI max height",
    "Fitted mz is outside the expected range",
    "Fitted mz is outside the expected range",
    "Fitted rt is outside the expected range",
    "Fitted rt is outside the expected range",
    "Fitted sigma_mz is too small",
    "Fitted sigma_rt is too small",
    "Fitted sigma_mz is too large",
    "Fitted sigma_rt is too large"
};

std::string Centroid::get_fit_failure_errors(const uint64_t &fit_failure_code, const std::vector<std::string> &error_messages) {
    // Function providing error messages for each bit in the fit_failure_code
    std::ostringstream result;
    result << "Fit Failure Errors:\n";

    // Iterate through each bit and check if it's set
    if (fit_failure_code == 0) {
        result << "No errors detected.\n";
    } else {
        for (size_t i = 0; i < error_messages.size(); ++i) {
            if (fit_failure_code & (1ULL << i)) {
                result << "- " << error_messages[i] << "\n";
            }
        }
    }

    return result.str();
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
    // log_itg("reset");

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
