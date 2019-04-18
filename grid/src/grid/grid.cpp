#include <cassert>
#include <cmath>

#include "grid/grid.hpp"
#include "utils/serialization.hpp"

// FIXME: Replace previous x_index.
uint64_t Grid::Mesh::x_index(double mz) {
    switch (this->instrument_type) {
        case Instrument::ORBITRAP: {
            double a = this->fwhm_mz / std::pow(this->reference_mz, 1.5);
            double b = (1 / std::sqrt(this->min_mz) - 1 / std::sqrt(mz));
            return static_cast<uint64_t>(this->k * 2 / a * b);
        } break;
        case Instrument::FTICR: {
            double a = 1 - this->min_mz / mz;
            double b = this->reference_mz * this->reference_mz;
            double c = this->fwhm_mz * this->min_mz;
            return static_cast<uint64_t>(this->k * a * b / c);
        } break;
        case Instrument::TOF: {
            return static_cast<uint64_t>(this->k * this->reference_mz /
                                         this->fwhm_mz *
                                         std::log(mz / this->min_mz));
        } break;
        case Instrument::QUAD: {
            // Same as the regular grid.
            return static_cast<uint64_t>(this->k * (mz - this->min_mz) /
                                         this->fwhm_mz);
        } break;
        case Instrument::UNKNOWN: {
            assert(false);  // Can't handle unknown instruments.
        } break;
    }
    assert(false);  // Can't handle unknown instruments.
    return 0;
}

// FIXME: Replace previous y_index.
uint64_t this_y_index(const RawData::RawData &raw_data, double rt, uint64_t k) {
    double delta_rt = raw_data.fwhm_rt / k;
    return std::ceil((rt - raw_data.min_rt) / delta_rt);
}

double Grid::Mesh::mz_at(uint64_t i) {
    switch (this->instrument_type) {
        case Instrument::ORBITRAP: {
            double a = 1 / std::sqrt(this->min_mz);
            double b = this->fwhm_mz / std::pow(this->reference_mz, 1.5) * i /
                       2 / this->k;
            double c = a - b;
            return 1 / (c * c);
        } break;
        case Instrument::FTICR: {
            double a = this->fwhm_mz * this->min_mz;
            double b = this->reference_mz * this->reference_mz;
            return this->min_mz / (1 - (a / b) * i / this->k);
        } break;
        case Instrument::TOF: {
            return this->min_mz *
                   std::exp(this->fwhm_mz / this->reference_mz * i / this->k);
        } break;
        case Instrument::QUAD: {
            double delta_mz = (this->max_mz - this->min_mz) /
                              static_cast<double>(this->n - 1);
            return this->min_mz + delta_mz * i / this->k;
        } break;
        case Instrument::UNKNOWN: {
            assert(false);  // Can't handle unknown instruments.
        } break;
    }
    assert(false);  // Can't handle unknown instruments.
    return 0;
}

double Grid::Mesh::rt_at(uint64_t j) {
    double delta_rt = (this->max_rt - this->min_rt) / (this->m - 1);
    return this->min_rt + delta_rt * j;
}

Grid::Mesh Grid::resample(const RawData::RawData &raw_data,
                          uint64_t num_samples_mz, uint64_t num_samples_rt,
                          double smoothing_coef_mz, double smoothing_coef_rt) {
    // Initialize the Mesh.
    Grid::Mesh mesh;
    mesh.k = num_samples_mz;
    mesh.t = num_samples_rt;
    mesh.reference_mz = raw_data.reference_mz;
    mesh.fwhm_mz = raw_data.reference_mz / raw_data.resolution_ms1;
    mesh.fwhm_rt = raw_data.fwhm_rt;
    mesh.instrument_type = raw_data.instrument_type;
    mesh.min_mz = raw_data.min_mz;
    mesh.max_mz = raw_data.max_mz;
    mesh.min_rt = raw_data.min_rt;
    mesh.max_rt = raw_data.max_rt;

    // Calculate the necessary dimensions for the Mesh.
    uint64_t n = mesh.x_index(raw_data.max_mz) + 1;
    uint64_t m = this_y_index(raw_data, raw_data.max_rt, num_samples_rt) + 1;
    mesh.n = n;
    mesh.m = m;
    mesh.data = std::vector<double>(n * m);
    mesh.bins_mz = std::vector<double>(n);
    mesh.bins_rt = std::vector<double>(m);

    // Generate bins_mz.
    for (size_t i = 0; i < n; ++i) {
        mesh.bins_mz[i] = mesh.mz_at(i);
    }

    // Generate bins_rt.
    for (size_t j = 0; j < m; ++j) {
        mesh.bins_rt[j] = mesh.rt_at(j);
    }

    // Pre-calculate the smoothing sigma values for all bins of the grid.
    double sigma_rt = RawData::fwhm_to_sigma(raw_data.fwhm_rt) *
                      smoothing_coef_rt / std::sqrt(2);
    auto sigma_mz_vec = std::vector<double>(n);
    for (size_t i = 0; i < n; ++i) {
        sigma_mz_vec[i] = RawData::fwhm_to_sigma(RawData::theoretical_fwhm(
                              raw_data, mesh.bins_mz[i])) *
                          smoothing_coef_mz / std::sqrt(2);
    }

    // Pre-calculate the kernel half widths for rt and all mzs.
    //
    // Since sigma_rt is constant, the size of the kernel will be the same
    // for the entire rt range.
    double delta_rt = (mesh.max_rt - mesh.min_rt) / (m - 1);
    uint64_t rt_kernel_hw = 3 * sigma_rt / delta_rt;
    auto mz_kernel_hw = std::vector<uint64_t>(n);
    for (size_t i = 0; i < n; ++i) {
        double sigma_mz = sigma_mz_vec[i];
        double delta_mz = 0;
        if (i == 0) {
            delta_mz = mesh.bins_mz[i + 1] - mesh.bins_mz[i];
        } else {
            delta_mz = mesh.bins_mz[i] - mesh.bins_mz[i - 1];
        }
        mz_kernel_hw[i] = 3 * sigma_mz / delta_mz;
    }

    // Gaussian splatting.
    {
        auto weights = std::vector<double>(n * m);
        for (size_t i = 0; i < raw_data.scans.size(); ++i) {
            const auto &scan = raw_data.scans[i];
            double current_rt = scan.retention_time;

            // Find the bin for the current retention time.
            size_t index_rt = this_y_index(raw_data, current_rt, mesh.t);

            // Find the min/max indexes for the rt kernel.
            size_t j_min = 0;
            if (index_rt >= rt_kernel_hw) {
                j_min = index_rt - rt_kernel_hw;
            }
            size_t j_max = mesh.m - 1;
            if ((index_rt + rt_kernel_hw) < mesh.m) {
                j_max = index_rt + rt_kernel_hw;
            }

            for (size_t k = 0; k < scan.num_points; ++k) {
                double current_intensity = scan.intensity[k];
                double current_mz = scan.mz[k];

                // Find the bin for the current mz.
                size_t index_mz = mesh.x_index(current_mz);

                double sigma_mz = sigma_mz_vec[index_mz];

                // Find the min/max indexes for the mz kernel.
                size_t i_min = 0;
                if (index_mz >= mz_kernel_hw[index_mz]) {
                    i_min = index_mz - mz_kernel_hw[index_mz];
                }
                size_t i_max = mesh.n - 1;
                if ((index_mz + mz_kernel_hw[index_mz]) < mesh.n) {
                    i_max = index_mz + mz_kernel_hw[index_mz];
                }

                for (size_t j = j_min; j <= j_max; ++j) {
                    for (size_t i = i_min; i <= i_max; ++i) {
                        double x = mesh.bins_mz[i];
                        double y = mesh.bins_rt[j];

                        // Calculate the Gaussian weight for this point.
                        double a = (x - current_mz) / sigma_mz;
                        double b = (y - current_rt) / sigma_rt;
                        double weight = std::exp(-0.5 * (a * a + b * b));

                        mesh.data[i + j * n] += weight * current_intensity;
                        weights[i + j * n] += weight;
                    }
                }
            }
        }
        for (size_t i = 0; i < (n * m); ++i) {
            double weight = weights[i];
            if (weight == 0) {
                weight = 1;
            }
            mesh.data[i] = mesh.data[i] / weight;
        }
    }

    // Gaussian smoothing.
    //
    // The Gaussian 2D filter is separable. We obtain the same result with
    // faster performance by applying two 1D kernel convolutions instead. This
    // is specially noticeable on the full image.
    {
        auto smoothed_data = std::vector<double>(mesh.n * mesh.m);

        // Retention time smoothing.
        for (size_t j = 0; j < mesh.m; ++j) {
            double current_rt = mesh.bins_rt[j];
            size_t min_k = 0;
            if (j >= rt_kernel_hw) {
                min_k = j - rt_kernel_hw;
            }
            size_t max_k = mesh.m - 1;
            if ((j + rt_kernel_hw) < mesh.m) {
                max_k = j + rt_kernel_hw;
            }
            for (size_t i = 0; i < mesh.n; ++i) {
                double sum_weights = 0;
                double sum_weighted_values = 0;
                for (size_t k = min_k; k <= max_k; ++k) {
                    double a = (current_rt - mesh.bins_rt[k]) / sigma_rt;
                    double weight = std::exp(-0.5 * (a * a));
                    sum_weights += weight;
                    sum_weighted_values += weight * mesh.data[i + k * mesh.n];
                }
                smoothed_data[i + j * mesh.n] =
                    sum_weighted_values / sum_weights;
            }
        }
        mesh.data = smoothed_data;
    }
    {
        auto smoothed_data = std::vector<double>(mesh.n * mesh.m);

        // mz smoothing.
        //
        // Since sigma_mz is not constant, we need to calculate the kernel half
        // width for each potential mz value.
        for (size_t i = 0; i < mesh.n; ++i) {
            double sigma_mz = sigma_mz_vec[i];
            double current_mz = mesh.bins_mz[i];

            size_t min_k = 0;
            if (i >= mz_kernel_hw[i]) {
                min_k = i - mz_kernel_hw[i];
            }
            size_t max_k = mesh.n - 1;
            if ((i + mz_kernel_hw[i]) < mesh.n) {
                max_k = i + mz_kernel_hw[i];
            }
            for (size_t j = 0; j < mesh.m; ++j) {
                double sum_weights = 0;
                double sum_weighted_values = 0;
                for (size_t k = min_k; k <= max_k; ++k) {
                    double a = (current_mz - mesh.bins_mz[k]) / sigma_mz;
                    double weight = std::exp(-0.5 * (a * a));
                    sum_weights += weight;
                    sum_weighted_values += weight * mesh.data[k + j * mesh.n];
                }
                smoothed_data[i + j * mesh.n] =
                    sum_weighted_values / sum_weights;
            }
        }
        mesh.data = smoothed_data;
    }
    return mesh;
}

// TODO(alex): Add optional normalization here.
bool Grid::splat(const Grid::Point &point, const Grid::Parameters &parameters,
                 std::vector<double> &data) {
    // For some instruments the peaks get wider in proportion to the adquired
    // mass. In order to maintain the same number of sampling points we will
    // scale the sigm_mz accordingly.
    double sigma_mz = Grid::sigma_mz(point.mz, parameters);
    double sigma_rt = Grid::sigma_rt(parameters);

    // Get the gaussian square dimensions. Note that we are using a square
    // kernel approximation for simplicity and computational efficiency.
    double min_rt = point.rt - 2 * sigma_rt;
    double max_rt = point.rt + 2 * sigma_rt;
    double min_mz = point.mz - 2 * sigma_mz;
    double max_mz = point.mz + 2 * sigma_mz;

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
            double x = Grid::mz_at(i, parameters);
            double y = Grid::rt_at(j, parameters);

            // Calculate the gaussian weight for this point.
            double a = (x - point.mz) / sigma_mz;
            double b = (y - point.rt) / sigma_rt;
            double weight = std::exp(-0.5 * (a * a + b * b));

            // Set the value, weight and counts.
            data[i + j * parameters.dimensions.n] += weight * point.value;
        }
    }
    return true;
}

double Grid::mz_at(uint64_t i, const Grid::Parameters &parameters) {
    // Warped grid.
    if (parameters.flags & Grid::Flags::WARPED_MESH) {
        switch (parameters.instrument_type) {
            case Instrument::ORBITRAP: {
                double a = 1 / std::sqrt(parameters.bounds.min_mz);
                double b = parameters.smoothing_params.sigma_mz /
                           std::pow(parameters.smoothing_params.mz, 1.5) * i;
                double c = (a - b);
                return 1 / (c * c);
            } break;
            case Instrument::FTICR: {
                double a = parameters.smoothing_params.sigma_mz *
                           parameters.bounds.min_mz;
                double b = parameters.smoothing_params.mz *
                           parameters.smoothing_params.mz;
                return parameters.bounds.min_mz / (1 - (a / b) * i);
            } break;
            case Instrument::TOF: {
                return parameters.bounds.min_mz *
                       std::exp(parameters.smoothing_params.sigma_mz /
                                parameters.smoothing_params.mz * i);
            } break;
            case Instrument::QUAD: {
                // Same as regular grid.
                double delta_mz =
                    (parameters.bounds.max_mz - parameters.bounds.min_mz) /
                    static_cast<double>(parameters.dimensions.n - 1);
                return parameters.bounds.min_mz + delta_mz * i;
            } break;
        }
    }

    // Regular grid.
    double delta_mz = (parameters.bounds.max_mz - parameters.bounds.min_mz) /
                      static_cast<double>(parameters.dimensions.n - 1);
    return parameters.bounds.min_mz + delta_mz * i;
}

double Grid::rt_at(uint64_t j, const Grid::Parameters &parameters) {
    double delta_rt = (parameters.bounds.max_rt - parameters.bounds.min_rt) /
                      static_cast<double>(parameters.dimensions.m - 1);
    return parameters.bounds.min_rt + delta_rt * j;
}

// TODO(alex): add unit tests for warped grid
uint64_t Grid::x_index(double mz, const Grid::Parameters &parameters) {
    // Regular grid.
    if (!(parameters.flags & Grid::Flags::WARPED_MESH)) {
        // In order to be consistent, the maximum value is mz + delta_mz. This
        // ensures all intervals contain the same number of points.
        double d = mz - parameters.bounds.min_mz;
        return static_cast<uint64_t>(d / parameters.smoothing_params.sigma_mz);
    }

    // Warped grid.
    switch (parameters.instrument_type) {
        case Instrument::ORBITRAP: {
            double a = std::pow(parameters.smoothing_params.mz, 1.5) /
                       parameters.smoothing_params.sigma_mz;
            double b = 1 / std::sqrt(parameters.bounds.min_mz);
            double c = std::sqrt(1 / mz);
            return static_cast<uint64_t>(a * (b - c));
        } break;
        case Instrument::FTICR: {
            double a = 1 - parameters.bounds.min_mz / mz;
            double b =
                parameters.smoothing_params.mz * parameters.smoothing_params.mz;
            double c =
                parameters.smoothing_params.sigma_mz * parameters.bounds.min_mz;
            return static_cast<uint64_t>(a * b / c);
        } break;
        case Instrument::TOF: {
            return static_cast<uint64_t>(
                parameters.smoothing_params.mz /
                parameters.smoothing_params.sigma_mz *
                std::log(mz / parameters.bounds.min_mz));
        } break;
        case Instrument::QUAD: {
            // Same as the regular grid.
            return static_cast<uint64_t>((mz - parameters.bounds.min_mz) /
                                         parameters.smoothing_params.sigma_mz);
        } break;
    }
    return 0;
}

uint64_t Grid::y_index(double rt, const Grid::Parameters &parameters) {
    // In order to be consistent, the maximum value is rt + delta_rt. This
    // ensures all intervals contain the same number of points.
    double d = rt - parameters.bounds.min_rt;
    return static_cast<uint64_t>(d / parameters.smoothing_params.sigma_rt);
}

double Grid::sigma_mz(double mz, const Grid::Parameters &parameters) {
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

double Grid::sigma_rt(const Grid::Parameters &parameters) {
    return parameters.smoothing_params.sigma_rt;
}

bool Grid::calculate_dimensions(Grid::Parameters &parameters) {
    if ((parameters.bounds.min_mz >= parameters.bounds.max_mz) ||
        (parameters.bounds.min_rt >= parameters.bounds.max_rt)) {
        return false;
    }

    auto n = Grid::x_index(parameters.bounds.max_mz, parameters);
    auto m = Grid::y_index(parameters.bounds.max_rt, parameters);
    parameters.dimensions = Grid::Dimensions{n + 1, m + 1};
    return true;
}
