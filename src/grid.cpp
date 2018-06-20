#include <cmath>

#include "grid.hpp"

bool Grid::splat(const Grid::Peak& peak, const Grid::Parameters& parameters,
                 std::vector<double>& data) {
    // For some instruments the peaks get wider in proportion to the adquired
    // mass. In order to maintain the same number of sampling points we will
    // scale the sigm_mz accordingly.
    double sigma_mz = Grid::sigma_mz(peak.mz, parameters);
    double sigma_rt = Grid::sigma_rt(parameters);

    // Get the gaussian square dimensions. Note that we are using a square
    // kernel approximation for simplicity and computational efficiency.
    double min_rt = peak.rt - 2 * sigma_rt;
    double max_rt = peak.rt + 2 * sigma_rt;
    double min_mz = peak.mz - 2 * sigma_mz;
    double max_mz = peak.mz + 2 * sigma_mz;

    if ((max_rt < parameters.bounds.min_rt) ||
        (max_mz < parameters.bounds.min_mz) ||
        (min_rt > parameters.bounds.max_rt) ||
        (min_mz > parameters.bounds.max_mz)) {
        return false;
    }

    // Even if the point lays outside the current grid, we still want to account
    // for it's contribution to the points in the frontier.
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

    for (size_t j = j_min; j <= j_max; ++j) {
        for (size_t i = i_min; i <= i_max; ++i) {
            // No need to do boundary check, since we are sure we are inside the
            // grid.
            double x = Grid::mz_at(i, parameters).value();
            double y = Grid::rt_at(j, parameters).value();

            // Calculate the gaussian weight for this point.
            double a = (x - peak.mz) / sigma_mz;
            double b = (y - peak.rt) / sigma_rt;
            double weight = std::exp(-0.5 * (a * a + b * b));

            // Set the value, weight and counts.
            auto old_value = Grid::value_at(i, j, parameters, data);
            Grid::set_value(i, j, old_value.value() + peak.value * weight,
                            parameters, data);
        }
    }
    return true;
}

std::optional<double> Grid::value_at(unsigned int i, unsigned int j,
                                     const Grid::Parameters& parameters,
                                     const std::vector<double>& data) {
    if (data.empty() || i > parameters.dimensions.n - 1 ||
        j > parameters.dimensions.m - 1) {
        return std::nullopt;
    }
    return data[i + j * parameters.dimensions.n];
}

bool Grid::set_value(unsigned int i, unsigned int j, double value,
                     const Grid::Parameters& parameters,
                     std::vector<double>& data) {
    if (data.empty() || i > parameters.dimensions.n - 1 ||
        j > parameters.dimensions.m - 1) {
        return false;
    }
    data[i + j * parameters.dimensions.n] = value;
    return true;
}

// TODO(alex): add unit tests for warped grid
std::optional<double> Grid::mz_at(unsigned int i,
                                  const Grid::Parameters& parameters) {
    if (parameters.dimensions.n * parameters.dimensions.m == 0 ||
        i > parameters.dimensions.n - 1) {
        return std::nullopt;
    }

    // Regular grid.
    if (!(parameters.flags & Grid::Flags::WARPED_MESH)) {
        auto delta_mz = (parameters.bounds.max_mz - parameters.bounds.min_mz) /
                        static_cast<double>(parameters.dimensions.n - 1);
        return parameters.bounds.min_mz + delta_mz * i;
    }

    // Warped grid.
    switch (parameters.instrument_type) {
        case Instrument::ORBITRAP: {
            double a = 1 / std::sqrt(parameters.bounds.min_mz);
            double b = parameters.smoothing_params.sigma_mz /
                       std::pow(parameters.smoothing_params.mz, 1.5) * i;
            double c = (a - b);
            return 1 / (c * c);
        } break;
        case Instrument::FTICR: {
            double a =
                parameters.smoothing_params.sigma_mz * parameters.bounds.min_mz;
            double b =
                parameters.smoothing_params.mz * parameters.smoothing_params.mz;
            return parameters.bounds.min_mz / (1 - (a / b) * i);
        } break;
        case Instrument::TOF: {
            return parameters.bounds.min_mz *
                   std::exp(parameters.smoothing_params.sigma_mz /
                            parameters.smoothing_params.mz * i);
        } break;
        case Instrument::QUAD: {
            // Same as regular grid.
            auto delta_mz =
                (parameters.bounds.max_mz - parameters.bounds.min_mz) /
                static_cast<double>(parameters.dimensions.n - 1);
            return parameters.bounds.min_mz + delta_mz * i;
        } break;
    }
    return std::nullopt;
}

std::optional<double> Grid::rt_at(unsigned int j,
                                  const Grid::Parameters& parameters) {
    if (parameters.dimensions.n * parameters.dimensions.m == 0 ||
        j > parameters.dimensions.m - 1) {
        return std::nullopt;
    }
    auto delta_rt = (parameters.bounds.max_rt - parameters.bounds.min_rt) /
                    static_cast<double>(parameters.dimensions.m - 1);
    return parameters.bounds.min_rt + delta_rt * j;
}

// TODO(alex): add unit tests for warped grid
unsigned int Grid::x_index(double mz, const Grid::Parameters& parameters) {
    // Regular grid.
    if (!(parameters.flags & Grid::Flags::WARPED_MESH)) {
        // In order to be consistent, the maximum value is mz + delta_mz. This
        // ensures all intervals contain the same number of points.
        double d = mz - parameters.bounds.min_mz;
        return static_cast<unsigned int>(d /
                                         parameters.smoothing_params.sigma_mz);
    }

    // Warped grid.
    switch (parameters.instrument_type) {
        case Instrument::ORBITRAP: {
            double a = std::pow(parameters.smoothing_params.mz, 1.5) /
                       parameters.smoothing_params.sigma_mz;
            double b = 1 / std::sqrt(parameters.bounds.min_mz);
            double c = std::sqrt(1 / mz);
            return static_cast<unsigned int>(a * (b - c));
        } break;
        case Instrument::FTICR: {
            double a = 1 - parameters.bounds.min_mz / mz;
            double b =
                parameters.smoothing_params.mz * parameters.smoothing_params.mz;
            double c =
                parameters.smoothing_params.sigma_mz * parameters.bounds.min_mz;
            return static_cast<unsigned int>(a * b / c);
        } break;
        case Instrument::TOF: {
            return static_cast<unsigned int>(
                parameters.smoothing_params.mz /
                parameters.smoothing_params.sigma_mz *
                std::log(mz / parameters.bounds.min_mz));
        } break;
        case Instrument::QUAD: {
            // Same as the regular grid.
            double d = mz - parameters.bounds.min_mz;
            return static_cast<unsigned int>(
                d / parameters.smoothing_params.sigma_mz);
        } break;
    }
    return 0;
}

unsigned int Grid::y_index(double rt, const Grid::Parameters& parameters) {
    // In order to be consistent, the maximum value is rt + delta_rt. This
    // ensures all intervals contain the same number of points.
    double d = rt - parameters.bounds.min_rt;
    return static_cast<unsigned int>(d / parameters.smoothing_params.sigma_rt);
}

double Grid::sigma_mz(double mz, const Grid::Parameters& parameters) {
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

double Grid::sigma_rt(const Grid::Parameters& parameters) {
    return parameters.smoothing_params.sigma_rt;
}

bool Grid::calculate_dimensions(Grid::Parameters& parameters) {
    if ((parameters.bounds.min_mz >= parameters.bounds.max_mz) ||
        (parameters.bounds.min_rt >= parameters.bounds.max_rt)) {
        return false;
    }

    auto n = Grid::x_index(parameters.bounds.max_mz, parameters);
    auto m = Grid::y_index(parameters.bounds.max_rt, parameters);
    parameters.dimensions = Grid::Dimensions{n + 1, m + 1};
    return true;
}
