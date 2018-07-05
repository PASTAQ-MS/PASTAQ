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

