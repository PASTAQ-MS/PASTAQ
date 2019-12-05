#ifndef CENTROID_CENTROID_HPP
#define CENTROID_CENTROID_HPP

#include <cstdint>
#include <optional>
#include <vector>

#include "grid/grid.hpp"

// In this namespace we can find the necessary functions to detect, quantify
// and manipulate peaks.
namespace Centroid {
// A LocalMax is defined as a point in the smoothed grid where the adjacent
// neighbours have a lower value than the current one.
struct LocalMax {
    double mz;
    double rt;
    double value;
};

// The detected Peak.
struct Peak {
    // ID of this peak. Should be kept for futher processing.
    uint64_t id;

    // Height, mz and rt values for the center of this peak (From the local
    // maxima coordinates on the grid).
    double local_max_mz;
    double local_max_rt;
    double local_max_height;

    // If the peak has been warped for retention time alignment, by adding the
    // following delta we can obtain the aligned retention time.
    double rt_delta;

    // Region of interest for this peak.
    double roi_min_mz;
    double roi_max_mz;
    double roi_min_rt;
    double roi_max_rt;

    // Simple estimation of the peak metrics on the raw data using the method
    // of moments.
    double raw_roi_mean_mz;
    double raw_roi_mean_rt;
    double raw_roi_sigma_mz;
    double raw_roi_sigma_rt;
    double raw_roi_skewness_mz;
    double raw_roi_skewness_rt;
    double raw_roi_kurtosis_mz;
    double raw_roi_kurtosis_rt;
    double raw_roi_max_height;
    double raw_roi_total_intensity;
    uint64_t raw_roi_num_points;
    uint64_t raw_roi_num_scans;

    // Gaussian fitting parameters.
    double fitted_height;
    double fitted_mz;
    double fitted_rt;
    double fitted_sigma_mz;
    double fitted_sigma_rt;
    double fitted_volume;
};

// Find all candidate points on the given grid by calculating the local maxima
// at each point of the grid. The local maxima is defined as follows: For the
// given indexes i and j the point at data[i][j] is greater than the neighbors
// in all 4 cardinal directions.
std::vector<LocalMax> find_local_maxima(const Grid::Grid &grid);

// Builds a Peak object for the given local_max.
std::optional<Peak> build_peak(const RawData::RawData &raw_data,
                               const LocalMax &local_max);

// Find the peaks in serial.
std::vector<Peak> find_peaks_serial(const RawData::RawData &raw_data,
                                    const Grid::Grid &grid, size_t max_peaks);

// Find the peaks in parallel.
std::vector<Peak> find_peaks_parallel(const RawData::RawData &raw_data,
                                      const Grid::Grid &grid, size_t max_peaks,
                                      size_t max_threads);

// Calculate the overlaping area between two peaks.
double peak_overlap(const Peak &peak_a, const Peak &peak_b);

// Calculate the cumulative similarity between two sets of peaks.
double cumulative_overlap(const std::vector<Peak> &set_a,
                          const std::vector<Peak> &set_b);

}  // namespace Centroid

#endif /* CENTROID_CENTROID_HPP */
