#include <algorithm>
#include <cmath>

#include "centroid.hpp"

std::vector<Centroid::Point> Centroid::find_local_maxima(
    const Grid::Parameters &parameters, const std::vector<double> &data) {
    std::vector<Centroid::Point> ret;
    uint64_t n_mz = parameters.dimensions.n;
    uint64_t n_rt = parameters.dimensions.m;
    for (size_t j = 1; j < n_rt; ++j) {
        for (size_t i = 1; i < n_mz; ++i) {
            int index = i + j * n_mz;

            // ---------------------------------------------------------
            // | top_left_value    | top_value    | top_right_value    |
            // ---------------------------------------------------------
            // | left_value        | value        | right_value        |
            // ---------------------------------------------------------
            // | bottom_left_value | bottom_value | bottom_right_value |
            // ---------------------------------------------------------
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
                peak.value = value;
                ret.push_back(peak);
            }
        }
    }
    return ret;
}

void Centroid::explore_peak_slope(uint64_t i, uint64_t j, double previous_value,
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
    // TODO(alex): Don't sort here or copy values, expect find_boundary to take
    // a list of already sorted peaks instead.
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

Centroid::Peak Centroid::build_peak(const Centroid::Point &local_max,
                                    const Grid::Parameters &parameters,
                                    const std::vector<double> &data) {
    Centroid::Peak peak = {};
    peak.i = local_max.i;
    peak.j = local_max.j;
    peak.height = local_max.value;
    peak.mz = Grid::mz_at(local_max.i, parameters);
    peak.rt = Grid::rt_at(local_max.j, parameters);

    // Extract the peak points and the boundary.
    Centroid::explore_peak_slope(local_max.i, local_max.j, -1, parameters, data,
                                 peak.points);

    // TODO(alex): Sort peaks here.
    // TODO(alex): error handling. What happens if the number of points is very
    // small? We should probably ignore peaks with less than 5 points so that it
    // has dimensionality in both mz and rt:
    //
    //   | |+| |
    //   |+|c|+|
    //   | |+| |
    peak.boundary = Centroid::find_boundary(peak.points);

    // Calculate the average background intensity from the boundary.
    {
        double boundary_sum = 0;
        for (const auto &point : peak.boundary) {
            boundary_sum += point.value;
        }
        peak.border_background = boundary_sum / peak.boundary.size();
    }

    // Calculate the total ion intensity on the peak for the values on the grid
    // and the sigma in mz and rt. The sigma is calculated by using the
    // algebraic formula for the variance of the random variable X:
    //
    //     Var(X) = E[X^2] - E[X]
    //
    // Where E[X] is the estimated value for X.
    //
    // In order to generalize this formula for the 2D blob, all values at the
    // same index will be aggregated together.
    //
    // TODO(alex): Note that this can cause catastrophic cancellation or
    // loss of significance. Probably the best option is to use a variant of
    // the Welford's method for computing the variance in a single pass. See:
    //
    //     http://jonisalonen.com/2013/deriving-welfords-method-for-computing-variance/
    //     https://ipfs.io/ipfs/QmXoypizjW3WknFiJnKLwHCnL72vedxjQkDDP1mXWo6uco/wiki/Algorithms_for_calculating_variance.html
    //
    {
        double height_sum = 0;
        double x_sum = 0;
        double y_sum = 0;
        double x_sig = 0;
        double y_sig = 0;
        for (const auto &point : peak.points) {
            double mz = Grid::mz_at(point.i, parameters);
            double rt = Grid::rt_at(point.j, parameters);

            height_sum += point.value;
            x_sum += point.value * mz;
            y_sum += point.value * rt;
            x_sig += point.value * mz * mz;
            y_sig += point.value * rt * rt;
        }
        peak.sigma_mz =
            std::sqrt((x_sig / height_sum) - std::pow(x_sum / height_sum, 2));
        peak.sigma_rt =
            std::sqrt((y_sig / height_sum) - std::pow(y_sum / height_sum, 2));

        // Make sure we don't have negative sigma values due to floating point
        // precision errors.
        peak.sigma_mz = peak.sigma_mz < 0 ? 1 : peak.sigma_mz;
        peak.sigma_rt = peak.sigma_rt < 0 ? 1 : peak.sigma_rt;

        peak.total_intensity = height_sum;
    }

    // Fit a weighted centroid to our data for the calculation of mz_centroid,
    // rt_centroid, height_centroid and total_intensity_centroid.
    // NOTE(alex): The original implementation uses +/- 4 points surrounding the
    // local maxima. Here we are using 2 * sigma_{rt,mz} to establish the size
    // of the window. We must verify that this does not have negative
    // implications on the performance and accuracy of the following tools in
    // the pipeline.
    {
        double sigma_rt = peak.sigma_rt;
        double sigma_mz = peak.sigma_mz;

        double min_rt = peak.rt - 2 * sigma_rt;
        double max_rt = peak.rt + 2 * sigma_rt;
        double min_mz = peak.mz - 2 * sigma_mz;
        double max_mz = peak.mz + 2 * sigma_mz;

        // Even if the point lays outside the current grid, we still want to
        // account for it's contribution to the points in the frontier.
        auto i_min = min_mz < parameters.bounds.min_mz
                         ? 0
                         : Grid::x_index(min_mz, parameters);
        auto j_min = min_rt < parameters.bounds.min_rt
                         ? 0
                         : Grid::y_index(min_rt, parameters);
        auto i_max = max_mz > parameters.bounds.max_mz
                         ? parameters.dimensions.n - 1
                         : Grid::x_index(max_mz, parameters);
        auto j_max = max_rt > parameters.bounds.max_rt
                         ? parameters.dimensions.m - 1
                         : Grid::y_index(max_rt, parameters);

        double height_sum = 0;
        double weights_sum = 0;
        double weighted_height_sum = 0;
        double x_sum = 0;
        double y_sum = 0;
        for (size_t j = j_min; j <= j_max; ++j) {
            for (size_t i = i_min; i <= i_max; ++i) {
                // No need to do boundary check, since we are sure we are inside
                // the grid.
                double mz = Grid::mz_at(i, parameters);
                double rt = Grid::rt_at(j, parameters);

                // Calculate the gaussian weight for this point.
                double a = (mz - peak.mz) / sigma_mz;
                double b = (rt - peak.rt) / sigma_rt;
                double weight = std::exp(-0.5 * (a * a + b * b));

                double height = data[i + j * parameters.dimensions.n];
                height_sum += height;
                weights_sum += weight;
                weighted_height_sum += height * weight;
                x_sum += height * mz;
                y_sum += height * rt;
            }
        }

        peak.mz_centroid = x_sum / height_sum;
        peak.rt_centroid = y_sum / height_sum;
        peak.total_intensity_centroid = height_sum;
        peak.height_centroid = weighted_height_sum / weights_sum;
    }

    return peak;
}
