#ifndef CENTROID_CENTROID_HPP
#define CENTROID_CENTROID_HPP

#include <cstdint>
#include <optional>
#include <tuple>
#include <vector>

#include "grid/grid.hpp"

namespace Centroid {

struct Parameters {
    uint64_t n_peaks;
    double threshold;
    Grid::Parameters grid_params;
};

// A Point represents the coordinates and value of those coordinates in the
// Grid.
struct Point {
    uint64_t i;
    uint64_t j;
    double value;
};

struct Peak {
    // ID of this peak. Should be kept for futher processing.
    size_t id;
    // Height,mz and rt values for the center of this peak (From the local
    // maxima coordinates on the mesh).
    double local_max_mz;
    double local_max_rt;
    double local_max_height;

    // Simple estimation of the peak metrics on the mesh values based on the
    // slope descent.
    //
    // Sumation of all intensities within the peak boundary. (Ignores holes,
    // i.e. does not interpolate values in case of non closed set).
    double slope_descent_total_intensity;
    // Estimated values for the position of the 2D peak based on the slope
    // descent points.
    double slope_descent_mean_mz;
    double slope_descent_mean_rt;
    // Estimated mz/rt values for the standard deviation of the peak in both
    // axes. (Ignores holes).
    double slope_descent_sigma_mz;
    double slope_descent_sigma_rt;
    // Average intensity on the boundary of the peak.
    double slope_descent_border_background;

    // Region of interest for this peak.
    double roi_min_mz;
    double roi_max_mz;
    double roi_min_rt;
    double roi_max_rt;
    // Simple estimation of the peak metrics on the raw data.
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
    uint64_t raw_roi_num_points_within_sigma;
    uint64_t raw_roi_num_scans;

    // FIXME: Better commments and documentation.
    // FIXME: method should take an enum instead of a string.
    // Extracted Ion Chromatogram. Can be performed with summation or maxima
    // (Total Ion Chromatogram/Base Peak Chromatogram). Returns two vectors
    // (retention_time, aggregated_intensity).
    std::tuple<std::vector<double>, std::vector<double>> xic(
        const RawData::RawData &raw_data, std::string method);
};

// Find all candidate points on the given mesh by calculating the local maxima
// at each point of the grid. The local maxima is defined as follows: For the
// given indexes i and j the point at data[i][j] is greater than the neighbors
// in all 8 directions. The resulting points are sorted in descending value
// order. If a maximum number of peaks is given on Centroid::Parameters.n_peaks,
// only the first n_peaks Points will be returned.
std::vector<Point> find_local_maxima(const Centroid::Parameters &parameters,
                                     const std::vector<double> &data);

// Find the boundary of the given bag of points.
std::vector<Point> find_boundary(std::vector<Point> &points);

// Find all points that belong to a given local max point via recursive local
// search of the slope of the peaks.
void explore_peak_slope(uint64_t i, uint64_t j, double previous_value,
                        const Centroid::Parameters &parameters,
                        const std::vector<double> &data,
                        std::vector<Point> &points);

// Builds a Peak object for the given local_max.
std::optional<Peak> build_peak(const Point &local_max,
                               const Centroid::Parameters &parameters,
                               const std::vector<double> &data);

}  // namespace Centroid

#endif /* CENTROID_CENTROID_HPP */
