#ifndef CENTROID_CENTROID_HPP
#define CENTROID_CENTROID_HPP

#include <cstdint>
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

// A Peak stores information regarding the peak detected in centroid. Note that
// for proper interpretation of the index values we need access to the
// Grid::Parameters from where these were calculated.
struct Peak {
    // Center of the peak in index space (Coordinates of local maxima).
    uint64_t i;
    uint64_t j;

    // Real mz/rt values for the center of this peak (From the local maxima
    // coordinates).
    double mz;
    double rt;
    // Height of the peak (Height of local maxima).
    double height;
    // Sumation of all intensities within the peak boundary. (Ignores holes,
    // i.e. does not interpolate values in case of non closed set).
    double total_intensity;

    // Estimated mz/rt values for the standard deviation of the peak in both
    // axes. (Ignores holes).
    double sigma_mz;
    double sigma_rt;

    // Metrics from the fitted weighted centroid.
    double mz_centroid;
    double rt_centroid;
    double height_centroid;
    double total_intensity_centroid;

    // Average intensity on the boundary of the peak.
    double border_background;

    // TODO(alex): Do we want to store these values?
    // Area under the curve of the extracted ion chromatogram. (Raw data)
    // double area_raw;
    // Area under the curve of the extracted ion chromatogram (Grid file,
    // smoothed data).
    // double area_smooth;
    // Sumation of all intensities within the peak boundary. (Raw data)
    // double total_ion_intensity_raw;
    // Number of ions contained within the peak boundary.
    // double number_measured_intensities;

    // All the points included within the peak boundary.
    std::vector<Point> points;
    std::vector<Point> boundary;
};

// Find all candidate points on the given mesh by calculating the local maxima
// at each point of the grid. The local maxima is defined as follows: For the
// given indexes i and j the point at data[i][j] is greater than the neighbors
// in all 8 directions. The resulting points are sorted in descending value
// order. If a maximum number of peaks is given on Centroid::Parameters.n_peaks,
// only the first n_peaks Points will be returned.
std::vector<Point> find_local_maxima(const Centroid::Parameters &parameters,
                                     const std::vector<double> &data);

// Find the boundary of the given bag of points. Note that this fuction sorts
// the points in place, so the original order is not preserved.
std::vector<Point> find_boundary(std::vector<Point> &points);

// Find all points that belong to a given local max point via recursive local
// search of the slope of the peaks.
void explore_peak_slope(uint64_t i, uint64_t j, double previous_value,
                        const Centroid::Parameters &parameters,
                        const std::vector<double> &data,
                        std::vector<Point> &points);

// Builds a Peak object for the given local_max.
Peak build_peak(const Point &local_max, const Centroid::Parameters &parameters,
                const std::vector<double> &data);

}  // namespace Centroid

#endif /* CENTROID_CENTROID_HPP */
