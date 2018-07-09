#include <algorithm>

#include "centroid.hpp"

std::vector<Centroid::Point> Centroid::find_local_maxima(
    const Grid::Parameters &parameters, const std::vector<double> &data) {
    std::vector<Centroid::Point> ret;
    unsigned int n_mz = parameters.dimensions.n;
    unsigned int n_rt = parameters.dimensions.m;
    for (size_t j = 1; j < n_rt; ++j) {
        for (size_t i = 1; i < n_mz; ++i) {
            int index = i + j * n_mz;

            double value = data[index];
            double right_value = data[index + 1];
            double left_value = data[index - 1];
            double top_value = data[index - n_mz];
            double top_right_value = data[index + 1 - n_mz];
            double top_left_value = data[index - 1 - n_mz];
            double bottom_value = data[index + n_mz];
            double bottom_right_value = data[index + 1 + n_mz];
            double bottom_left_value = data[index - 1 + n_mz];

            if ((value != 0) && (value > left_value) && (value > right_value) &&
                (value > top_value) && (value > top_right_value) &&
                (value > top_left_value) && (value > bottom_value) &&
                (value > bottom_left_value) && (value > bottom_right_value)) {
                Centroid::Point peak = {};
                peak.i = i;
                peak.j = j;
                peak.height = value;
                ret.push_back(peak);
            }
        }
    }
    return ret;
}

// FIXME(alex): This does not account for the case where a peak might be
// contained in another peak. Local search is necessary in that case.
std::vector<Centroid::Point> Centroid::find_peak_points(
    const Centroid::Point &point, const Grid::Parameters &parameters,
    const std::vector<double> &data) {
    // Initialize return vector.
    std::vector<Centroid::Point> ret;
    ret.push_back(point);

    // TODO(alex): Set to a greater value, select by user or calculate it (For
    // example accounting for maximum floating point precision.
    double threshold = 0;

    unsigned int max_i = 0;
    unsigned int min_i = 0;
    unsigned int max_j = 0;
    unsigned int min_j = 0;

    // Find right boundary.
    {
        double previous_value = point.height;
        for (unsigned int i = point.i; i < parameters.dimensions.n; ++i) {
            int index = i + point.j * parameters.dimensions.n;
            double value = data[index];
            if (value <= threshold || value > previous_value) {
                max_i = i;
                break;
            }
            if (i != point.i) {
                ret.push_back({i, point.j, value});
            }
            previous_value = value;
        }
    }
    // Find left boundary.
    {
        double previous_value = point.height;
        for (unsigned int i = point.i; i + 1 > 0; --i) {
            int index = i + point.j * parameters.dimensions.n;
            double value = data[index];
            if (value <= threshold || value > previous_value) {
                min_i = i;
                break;
            }
            if (i != point.i) {
                ret.push_back({i, point.j, value});
            }
            previous_value = value;
        }
    }
    // Find bottom boundary.
    {
        double previous_value = point.height;
        for (unsigned int j = point.j; j < parameters.dimensions.m; ++j) {
            int index = point.i + j * parameters.dimensions.n;
            double value = data[index];
            if (value <= threshold || value > previous_value) {
                max_j = j;
                break;
            }
            if (j != point.j) {
                ret.push_back({point.i, j, value});
            }
            previous_value = value;
        }
    }
    // Find top boundary.
    {
        double previous_value = point.height;
        for (unsigned int j = point.j; j + 1 > 0; --j) {
            int index = point.i + j * parameters.dimensions.n;
            double value = data[index];
            if (value <= threshold || value > previous_value) {
                min_j = j;
                break;
            }
            if (j != point.j) {
                ret.push_back({point.i, j, value});
            }
            previous_value = value;
        }
    }

    // Peak quadrants:
    //
    // |----|----|
    // | Q4 | Q1 |
    // |----|----|
    // | Q3 | Q2 |
    // |----|----|
    //
    // Find Q1 points
    for (unsigned int j = point.j - 1; j + 1 > min_j; --j) {
        double previous_value = data[point.i + j * parameters.dimensions.n];
        for (unsigned int i = point.i + 1; i < max_i; ++i) {
            int index = i + j * parameters.dimensions.n;
            double value = data[index];
            if (value <= threshold || value > previous_value) {
                break;
            }
            if (i != point.i || j != point.j) {
                ret.push_back({i, j, value});
            }
            previous_value = value;
        }
    }
    // Find Q2 points
    for (unsigned int j = point.j + 1; j < max_j; ++j) {
        double previous_value = data[point.i + j * parameters.dimensions.n];
        for (unsigned int i = point.i + 1; i < max_i; ++i) {
            int index = i + j * parameters.dimensions.n;
            double value = data[index];
            if (value <= threshold || value > previous_value) {
                break;
            }
            if (i != point.i || j != point.j) {
                ret.push_back({i, j, value});
            }
            previous_value = value;
        }
    }
    // Find Q3 points
    for (unsigned int j = point.j - 1; j + 1 > min_j; --j) {
        double previous_value = data[point.i + j * parameters.dimensions.n];
        for (unsigned int i = point.i - 1; i + 1 > min_i; --i) {
            int index = i + j * parameters.dimensions.n;
            double value = data[index];
            if (value <= threshold || value > previous_value) {
                break;
            }
            if (i != point.i || j != point.j) {
                ret.push_back({i, j, value});
            }
            previous_value = value;
        }
    }
    // Find Q4 points
    for (unsigned int j = point.j + 1; j < max_j; ++j) {
        double previous_value = data[point.i + j * parameters.dimensions.n];
        for (unsigned int i = point.i - 1; i + 1 > min_i; --i) {
            int index = i + j * parameters.dimensions.n;
            double value = data[index];
            if (value <= threshold || value > previous_value) {
                break;
            }
            if (i != point.i || j != point.j) {
                ret.push_back({i, j, value});
            }
            previous_value = value;
        }
    }

    return ret;
}

std::vector<Centroid::Point> Centroid::find_boundary(
    const std::vector<Centroid::Point> &points) {
    size_t n_points = points.size();
    // Creates a copy of the input vector to avoid modifying the original.
    std::vector<Centroid::Point> points_copy(n_points);
    std::copy(points.begin(), points.end(), points_copy.begin());

    // Under the constraints of the grid coordinates, we need at least 5 points
    // in order to have a boundary that does not contain all the points in the
    // initial set.
    if (n_points < 5) {
        return points_copy;
    }

    // Sort the points by y and then x coordinates.
    auto sort_points = [](const Centroid::Point &p1,
                          const Centroid::Point &p2) -> bool {
        return (p1.j < p2.j) || ((p1.j == p2.j) && (p1.i < p2.i));
    };
    std::stable_sort(points_copy.begin(), points_copy.end(), sort_points);

    // Peeking the first and last points to get the minimum and maximum row
    // number.
    size_t min_j = points_copy[0].j;
    size_t max_j = points_copy[n_points - 1].j;

    // Extract the boundary points through row major iteration of the virtual
    // grid created by the bag of points.
    std::vector<Centroid::Point> boundary_points;
    for (size_t k = 0; k < n_points; ++k) {
        auto point = points_copy[k];
        if (point.j == min_j || point.j == max_j) {
            boundary_points.push_back(point);
        } else {
            boundary_points.push_back(point);
            // Find last element in the row.
            while (k + 1 < n_points && points_copy[k + 1].j == point.j) {
                ++k;
            }
            boundary_points.push_back(points_copy[k]);
        }
    }
    return boundary_points;
}

void Centroid::explore_peak_slope(unsigned int i, unsigned int j,
                                  double previous_value,
                                  const Grid::Parameters &parameters,
                                  const std::vector<double> &data,
                                  std::vector<Centroid::Point> &points) {
    // Check that the point has not being already included.
    for (const auto &point : points) {
        if (point.i == i && point.j == j) {
            return;
        }
    }

    // here set threshold to MAXIMUM of fractional peak height and a min
    // value
    // double bestthresh = thresh * pheight;
    // if (bestthresh < peakheightmin) bestthresh = peakheightmin;
    // TODO(alex): Set to a greater value, select by user or calculate it (For
    // example accounting for maximum floating point precision.
    double threshold = 0;

    double value = data[i + j * parameters.dimensions.n];
    if (previous_value >= 0 && (previous_value < value || value <= threshold)) {
        return;
    }

    points.push_back({i, j, value});

    // Return if we are at the edge of the grid.
    // FIXME(alex): Can this cause problems if we could continue exploring
    // downwards?
    if (i < 1 || i >= parameters.dimensions.n - 1 || j < 1 ||
        j >= parameters.dimensions.m - 1) {
        return;
    }

    Centroid::explore_peak_slope(i - 1, j, value, parameters, data, points);
    Centroid::explore_peak_slope(i + 1, j, value, parameters, data, points);
    Centroid::explore_peak_slope(i, j + 1, value, parameters, data, points);
    Centroid::explore_peak_slope(i, j - 1, value, parameters, data, points);
    Centroid::explore_peak_slope(i - 1, j - 1, value, parameters, data, points);
    Centroid::explore_peak_slope(i + 1, j + 1, value, parameters, data, points);
    Centroid::explore_peak_slope(i - 1, j + 1, value, parameters, data, points);
    Centroid::explore_peak_slope(i + 1, j - 1, value, parameters, data, points);
}
