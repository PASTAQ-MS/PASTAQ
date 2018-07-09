#ifndef CENTROID_CENTROID_HPP
#define CENTROID_CENTROID_HPP

#include <vector>

#include "grid.hpp"

namespace Centroid {

struct Point {
    unsigned int i;
    unsigned int j;
    double height;
};

struct Peak {
    // Center of the peak in index space.
    unsigned int i;
    unsigned int j;

    // Real mz/rt values for the center of this peak.
    double mz;
    double rt;

    // Real mz/rt values for the standard deviation of the peak in both axis.
    double sigma_mz;
    double sigma_rt;

    // TODO(alex): Do we want to store these values? If so, should we calculate
    // this over the height or the height/noise?
    // double fwhm_mz;
    // double fwhm_rt;

    // Height of the peak centroid.
    double height;
    // Sumation of all intensities within the peak boundary.
    double total_intensity;
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

// Find all points that belong to a given local max point.
std::vector<Point> find_peak_points(const Point &point,
                                    const Grid::Parameters &parameters,
                                    const std::vector<double> &data);

// Find the boundary of the given bag of points.
std::vector<Point> find_boundary(const std::vector<Point> &points);

void explore_peak_slope(unsigned int i, unsigned int j, double previous_value,
                        const Grid::Parameters &parameters,
                        const std::vector<double> &data,
                        std::vector<Point> &points);

Peak build_peak(const Point &local_max, const Grid::Parameters &parameters,
                const std::vector<double> &data);
}  // namespace Centroid

#endif /* CENTROID_CENTROID_HPP */
