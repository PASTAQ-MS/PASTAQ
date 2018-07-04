#ifndef CENTROID_CENTROID_HPP
#define CENTROID_CENTROID_HPP

#include <vector>

#include "grid.hpp"

namespace Centroid {

struct Peak {
    // TODO(alex): Do we want to store both index AND values? If not, which ones
    // do we keep?
    unsigned int i;
    unsigned int j;
    // TODO(alex): Do we want to store this or calculate them using
    // Grid::Parameters?
    // double mz;
    // double rt;
    // TODO(alex): How many more of these do we want? VCentroid? XSigma? YSigma?
    // Point count? S/N ratio? Background noise?
    double height;
    // double volume;
};

struct Point {
    unsigned int i;
    unsigned int j;
    double height;
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

}  // namespace Centroid

#endif /* CENTROID_CENTROID_HPP */
