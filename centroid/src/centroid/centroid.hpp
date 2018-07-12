#ifndef CENTROID_CENTROID_HPP
#define CENTROID_CENTROID_HPP

#include <vector>

#include "grid/grid.hpp"

// TODO(alex): Add doc.
namespace Centroid {

// TODO(alex): Add doc.
struct Point {
    unsigned int i;
    unsigned int j;
    double height;
};

// TODO(alex): Add doc.
struct Peak {
    // Center of the peak in index space (Coordinates of local maxima).
    unsigned int i;
    unsigned int j;

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

// Find all points that belong to a given local max point.
std::vector<Point> find_peak_points(const Point &point,
                                    const Grid::Parameters &parameters,
                                    const std::vector<double> &data);

// Find the boundary of the given bag of points.
std::vector<Point> find_boundary(const std::vector<Point> &points);

// TODO(alex): Add doc.
void explore_peak_slope(unsigned int i, unsigned int j, double previous_value,
                        const Grid::Parameters &parameters,
                        const std::vector<double> &data,
                        std::vector<Point> &points);

// TODO(alex): Add doc.
Peak build_peak(const Point &local_max, const Grid::Parameters &parameters,
                const std::vector<double> &data);

}  // namespace Centroid

#endif /* CENTROID_CENTROID_HPP */
