#ifndef CENTROID_CENTROID_HPP
#define CENTROID_CENTROID_HPP

#include <cstdint>
#include <vector>

#include "grid/grid.hpp"

namespace Centroid {

// A Point represents the coordinates and value of those coordinates in the
// Grid.
struct Point {
    // TODO(alex): should we store instead the mz/rt directly to remove the
    // Grid::Parameters dependency?
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

    // TODO(alex): Do we want to store these values? If so, should we calculate
    // this over the height or the height/noise?
    // double fwhm_mz;
    // double fwhm_rt;
    // TODO(alex): Store SN values? sn_height or sn_centroid

    // Average intensity on the boundary of the peak.
    double border_background;

    // TODO(alex): Do we want to store these values?
    // Volume of the 3D peak.
    // double volume;
    // Area under the curve of the extracted ion chromatogram. (Raw file)
    // double area_raw;
    // Area under the curve of the extracted ion chromatogram (Grid file,
    // smoothed data).
    // double area_smooth;
    // Sumation of all intensities within the peak boundary. (Grid file,
    // smoothed)
    // double total_ion_intensity_smooth;
    // Number of ions contained within the peak boundary.
    // double number_measured_intensities;

    // All the points included within the peak boundary.
    std::vector<Point> points;
    std::vector<Point> boundary;
};

// Find all candidate points on the given mesh by calculating the local maxima
// at each point of the grid. The local maxima is defined as follows: For the
// given indexes i and j the point at data[i][j] is greater than the neighbors
// in all 8 directions.
std::vector<Point> find_local_maxima(const Grid::Parameters &parameters,
                                     const std::vector<double> &data);

// Find the boundary of the given bag of points.
std::vector<Point> find_boundary(const std::vector<Point> &points);

// Find all points that belong to a given local max point via recursive local
// search of the slope of the peaks.
void explore_peak_slope(uint64_t i, uint64_t j, double previous_value,
                        const Grid::Parameters &parameters,
                        const std::vector<double> &data,
                        std::vector<Point> &points);

// Builds a Peak object for the given local_max.
Peak build_peak(const Point &local_max, const Grid::Parameters &parameters,
                const std::vector<double> &data);

}  // namespace Centroid

#endif /* CENTROID_CENTROID_HPP */
