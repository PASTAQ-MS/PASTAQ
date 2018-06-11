#include <cmath>

#include "tapp-grid.hpp"

bool Grid::splat(double mz, double rt, double value,
                 Grid::Parameters& parameters, std::vector<double>& data) {
    // For some instruments the peaks get wider in proportion to the adquired
    // mass. In order to maintain the same number of sampling points we will
    // scale the sigm_mz accordingly.
    double sigma_mz = Grid::sigma_mz(mz, parameters);
    double sigma_rt = Grid::sigma_rt(parameters);

    // Get the gaussian square dimensions. Note that we are using a square
    // kernel approximation for simplicity and computational efficiency.
    double min_rt = rt - 2 * sigma_rt;
    double max_rt = rt + 2 * sigma_rt;
    double min_mz = mz - 2 * sigma_mz;
    double max_mz = mz + 2 * sigma_mz;

    if ((max_rt < parameters.bounds.min_rt) ||
        (max_mz < parameters.bounds.min_mz) ||
        (min_rt > parameters.bounds.max_rt) ||
        (min_mz > parameters.bounds.max_mz)) {
        return false;
    }

    // Even if the point lays outside the current grid, we still want to account
    // for it's contribution to the points in the frontier.
    auto i_min = Grid::x_index(min_mz, parameters)
                     ? Grid::x_index(min_mz, parameters).value()
                     : 0;
    auto j_min = Grid::y_index(min_rt, parameters)
                     ? Grid::y_index(min_rt, parameters).value()
                     : 0;
    auto i_max = Grid::x_index(max_mz, parameters)
                     ? Grid::x_index(max_mz, parameters).value()
                     : parameters.dimensions.n - 1;
    auto j_max = Grid::y_index(max_rt, parameters)
                     ? Grid::y_index(max_rt, parameters).value()
                     : parameters.dimensions.m - 1;

    double x0 = mz;
    double y0 = rt;
    for (size_t j = j_min; j <= j_max; ++j) {
        for (size_t i = i_min; i <= i_max; ++i) {
            // No need to do boundary check, since we are sure we are inside the
            // grid.
            double x = Grid::mz_at(i, parameters).value();
            double y = Grid::rt_at(j, parameters).value();

            // Calculate the gaussian weight for this point.
            double a = (x - x0) / sigma_mz;
            double b = (y - y0) / sigma_rt;
            double weight = std::exp(-0.5 * (a * a + b * b));

            // Set the value, weight and counts.
            auto old_value = Grid::value_at(i, j, parameters, data);
            Grid::set_value(i, j, old_value.value() + value * weight,
                            parameters, data);
        }
    }
    return true;
}

std::optional<double> Grid::value_at(unsigned int i, unsigned int j,
                                     Grid::Parameters& parameters,
                                     std::vector<double> data) {
    if (data.empty() || i > parameters.dimensions.n - 1 ||
        j > parameters.dimensions.m - 1) {
        return std::nullopt;
    }
    return data[i + j * parameters.dimensions.n];
}

bool Grid::set_value(unsigned int i, unsigned int j, double value,
                     Grid::Parameters& parameters, std::vector<double>& data) {
    if (data.empty() || i > parameters.dimensions.n - 1 ||
        j > parameters.dimensions.m - 1) {
        return false;
    }
    data[i + j * parameters.dimensions.n] = value;
    return true;
}

std::optional<double> Grid::mz_at(unsigned int i,
                                  Grid::Parameters& parameters) {
    if (parameters.dimensions.n * parameters.dimensions.m == 0 ||
        i > parameters.dimensions.n - 1) {
        return std::nullopt;
    }

    auto delta_mz = (parameters.bounds.max_mz - parameters.bounds.min_mz) /
                    static_cast<double>(parameters.dimensions.n - 1);
    return parameters.bounds.min_mz + delta_mz * i;
}

std::optional<double> Grid::rt_at(unsigned int j,
                                  Grid::Parameters& parameters) {
    if (parameters.dimensions.n * parameters.dimensions.m == 0 ||
        j > parameters.dimensions.m - 1) {
        return std::nullopt;
    }
    auto delta_rt = (parameters.bounds.max_rt - parameters.bounds.min_rt) /
                    static_cast<double>(parameters.dimensions.m - 1);
    return parameters.bounds.min_rt + delta_rt * j;
}

std::optional<unsigned int> Grid::x_index(double mz,
                                          Grid::Parameters& parameters) {
    // In order to be consistent, the maximum value is mz + delta_mz. This
    // ensures all intervals contain the same number of points.
    double delta_mz = (parameters.bounds.max_mz - parameters.bounds.min_mz) /
                      static_cast<double>(parameters.dimensions.n - 1);
    if (mz < parameters.bounds.min_mz ||
        mz > parameters.bounds.max_mz + delta_mz) {
        return std::nullopt;
    }
    double d = mz - parameters.bounds.min_mz;
    auto i = static_cast<unsigned int>(d / delta_mz);
    if (i > parameters.dimensions.n - 1) {
        return std::nullopt;
    }
    return i;
}

std::optional<unsigned int> Grid::y_index(double rt,
                                          Grid::Parameters& parameters) {
    // In order to be consistent, the maximum value is rt + delta_rt. This
    // ensures all intervals contain the same number of points.
    double delta_rt = (parameters.bounds.max_rt - parameters.bounds.min_rt) /
                      static_cast<double>(parameters.dimensions.m - 1);
    if (rt < parameters.bounds.min_rt ||
        rt > parameters.bounds.max_rt + delta_rt) {
        return std::nullopt;
    }
    double d = rt - parameters.bounds.min_rt;
    auto j = static_cast<unsigned int>(d / delta_rt);
    if (j > parameters.dimensions.m - 1) {
        return std::nullopt;
    }
    return j;
}

double Grid::sigma_mz(double mz, Grid::Parameters& parameters) {
    double sigma_mz = 0.0;
    switch (parameters.instrument_type) {
        case Instrument::ORBITRAP: {
            sigma_mz = parameters.smoothing_params.sigma_mz *
                       std::pow(mz / parameters.smoothing_params.mz, 1.5);
        } break;
        case Instrument::FTICR: {
            sigma_mz = parameters.smoothing_params.sigma_mz *
                       std::pow(mz / parameters.smoothing_params.mz, 2);
        } break;
        case Instrument::TOF: {
            sigma_mz = parameters.smoothing_params.sigma_mz * mz /
                       parameters.smoothing_params.mz;
        } break;
        case Instrument::QUAD: {
            // QUAD/IONTRAP instruments maintain the same resolution across all
            // mass range.
            sigma_mz = parameters.smoothing_params.sigma_mz;
        } break;
    }
    return sigma_mz;
}

double Grid::sigma_rt(Grid::Parameters& parameters) {
    return parameters.smoothing_params.sigma_rt;
}
