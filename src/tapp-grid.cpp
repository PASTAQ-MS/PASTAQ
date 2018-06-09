#include <cmath>

#include "tapp-grid.hpp"

bool Grid::splash(Grid::Interface &grid, double mz, double rt, double value,
                  double sigma_rt) {
    // For some instruments the peaks get wider in proportion to the adquired
    // mass. In order to maintain the same number of sampling points we will
    // scale the sigm_mz accordingly.
    double sigma_mz = grid.sigma_at_mz(mz);

    // Get the gaussian square dimensions. Note that we are using a square
    // kernel approximation for simplicity and computational efficiency.
    double min_rt = rt - 2 * sigma_rt;
    double max_rt = rt + 2 * sigma_rt;
    double min_mz = mz - 2 * sigma_mz;
    double max_mz = mz + 2 * sigma_mz;

    if ((max_rt < grid.bounds().min_rt) || (max_mz < grid.bounds().min_mz) ||
        (min_rt > grid.bounds().max_rt) || (min_mz > grid.bounds().max_mz)) {
        return false;
    }

    // Even if the point lays outside the current grid, we still want to account
    // for it's contribution to the points in the frontier.
    auto i_min = grid.x_index(min_mz) ? grid.x_index(min_mz).value() : 0;
    auto j_min = grid.y_index(min_rt) ? grid.y_index(min_rt).value() : 0;
    auto i_max =
        grid.x_index(max_mz) ? grid.x_index(max_mz).value() : grid.dim().n - 1;
    auto j_max =
        grid.y_index(max_rt) ? grid.y_index(max_rt).value() : grid.dim().m - 1;

    double x0 = mz;
    double y0 = rt;
    for (unsigned int j = j_min; j <= j_max; ++j) {
        for (unsigned int i = i_min; i <= i_max; ++i) {
            // No need to do boundary check, since we are sure we are inside the
            // grid.
            double x = grid.mz_at(i).value();
            double y = grid.rt_at(j).value();

            // Calculate the gaussian weight for this point.
            double a = (x - x0) / sigma_mz;
            double b = (y - y0) / sigma_rt;
            double weight = std::exp(-0.5 * (a * a + b * b));

            // Set the value, weight and counts.
            auto old_value = grid.value_at(i, j);
            grid.set_value(i, j, old_value.value() + value * weight);
        }
    }
    return true;
}
