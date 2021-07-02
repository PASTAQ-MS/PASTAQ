#include <cassert>
#include <cmath>

#include "grid/grid.hpp"
#include "utils/serialization.hpp"
#include "utils/search.hpp"

uint64_t Grid::x_index(const Grid &grid, double mz) {
    switch (grid.instrument_type) {
        case Instrument::ORBITRAP: {
            double a = grid.fwhm_mz / std::pow(grid.reference_mz, 1.5);
            double b = (1 / std::sqrt(grid.min_mz) - 1 / std::sqrt(mz));
            return static_cast<uint64_t>(grid.k * 2 / a * b);
        } break;
        case Instrument::FTICR: {
            double a = 1 - grid.min_mz / mz;
            double b = grid.reference_mz * grid.reference_mz;
            double c = grid.fwhm_mz * grid.min_mz;
            return static_cast<uint64_t>(grid.k * a * b / c);
        } break;
        case Instrument::TOF: {
            return static_cast<uint64_t>(grid.k * grid.reference_mz /
                                         grid.fwhm_mz *
                                         std::log(mz / grid.min_mz));
        } break;
        case Instrument::QUAD: {
            // Same as the regular grid.
            return static_cast<uint64_t>(grid.k * (mz - grid.min_mz) /
                                         grid.fwhm_mz);
        } break;
        default: {
            assert(false);  // Can't handle unknown instruments.
            return 0;
        } break;
    }
}

uint64_t Grid::y_index(const Grid &grid, double rt) {
    double delta_rt = grid.fwhm_rt / grid.k;
    return std::ceil((rt - grid.min_rt) / delta_rt);
}

double Grid::mz_at(const Grid &grid, uint64_t i) {
    switch (grid.instrument_type) {
        case Instrument::ORBITRAP: {
            double a = 1 / std::sqrt(grid.min_mz);
            double b = grid.fwhm_mz / std::pow(grid.reference_mz, 1.5) * i / 2 /
                       grid.k;
            double c = a - b;
            return 1 / (c * c);
        } break;
        case Instrument::FTICR: {
            double a = grid.fwhm_mz * grid.min_mz;
            double b = grid.reference_mz * grid.reference_mz;
            return grid.min_mz / (1 - (a / b) * i / grid.k);
        } break;
        case Instrument::TOF: {
            return grid.min_mz *
                   std::exp(grid.fwhm_mz / grid.reference_mz * i / grid.k);
        } break;
        case Instrument::QUAD: {
            double delta_mz =
                (grid.max_mz - grid.min_mz) / static_cast<double>(grid.n - 1);
            return grid.min_mz + delta_mz * i / grid.k;
        } break;
        default: {
            assert(false);  // Can't handle unknown instruments.
            return 0;
        } break;
    }
}

double Grid::rt_at(const Grid &grid, uint64_t j) {
    double delta_rt = (grid.max_rt - grid.min_rt) / (grid.m - 1);
    return grid.min_rt + delta_rt * j;
}

Grid::Grid Grid::resample(const RawData::RawData &raw_data,
                          const ResampleParams &params) {
    // Initialize the Grid.
    Grid grid;
    grid.k = params.num_samples_mz;
    grid.t = params.num_samples_rt;
    grid.reference_mz = raw_data.reference_mz;
    grid.fwhm_mz = raw_data.reference_mz / raw_data.resolution_ms1;
    grid.fwhm_rt = raw_data.fwhm_rt;
    grid.instrument_type = raw_data.instrument_type;
    grid.min_mz = raw_data.min_mz;
    grid.max_mz = raw_data.max_mz;
    grid.min_rt = raw_data.min_rt;
    grid.max_rt = raw_data.max_rt;

    // Calculate the necessary dimensions for the Grid.
    uint64_t n = x_index(grid, raw_data.max_mz) + 1;
    uint64_t m = y_index(grid, raw_data.max_rt) + 1;
    grid.n = n;
    grid.m = m;
    grid.data = std::vector<double>(n * m);
    grid.bins_mz = std::vector<double>(n);
    grid.bins_rt = std::vector<double>(m);

    // Generate bins_mz.
    for (size_t i = 0; i < n; ++i) {
        grid.bins_mz[i] = mz_at(grid, i);
    }

    // Generate bins_rt.
    for (size_t j = 0; j < m; ++j) {
        grid.bins_rt[j] = rt_at(grid, j);
    }

    // Pre-calculate the smoothing sigma values for all bins of the grid.
    double sigma_rt = RawData::fwhm_to_sigma(raw_data.fwhm_rt) *
                      params.smoothing_coef_rt / std::sqrt(2);
    auto sigma_mz_vec = std::vector<double>(n);
    for (size_t i = 0; i < n; ++i) {
        sigma_mz_vec[i] = RawData::fwhm_to_sigma(RawData::theoretical_fwhm(
                              raw_data, grid.bins_mz[i])) *
                          params.smoothing_coef_mz / std::sqrt(2);
    }

    // Pre-calculate the kernel half widths for rt and mz.
    //
    // Since sigma_rt is constant, the size of the kernel will be the same
    // for the entire rt range. The same applies for mz, as we are using a
    // warped grid that keeps the number of sampling points constant across the
    // entire m/z range.
    double delta_rt = grid.fwhm_rt / params.num_samples_rt;
    double delta_mz = grid.fwhm_mz / params.num_samples_mz;
    double sigma_mz_ref = RawData::fwhm_to_sigma(grid.fwhm_mz);
    uint64_t rt_kernel_hw = 3 * sigma_rt / delta_rt;
    uint64_t mz_kernel_hw = 3 * sigma_mz_ref / delta_mz;

    // Gaussian splatting.
    {
        auto weights = std::vector<double>(n * m);
        for (size_t s = 0; s < raw_data.scans.size(); ++s) {
            const auto &scan = raw_data.scans[s];
            double current_rt = scan.retention_time;

            // Find the bin for the current retention time.
            size_t index_rt = y_index(grid, current_rt);

            // Find the min/max indexes for the rt kernel.
            size_t j_min = 0;
            if (index_rt >= rt_kernel_hw) {
                j_min = index_rt - rt_kernel_hw;
            }
            size_t j_max = grid.m - 1;
            if ((index_rt + rt_kernel_hw) < grid.m) {
                j_max = index_rt + rt_kernel_hw;
            }

            for (size_t k = 0; k < scan.num_points; ++k) {
                double current_intensity = scan.intensity[k];
                double current_mz = scan.mz[k];

                // Find the bin for the current mz.
                size_t index_mz = x_index(grid, current_mz);

                double sigma_mz = sigma_mz_vec[index_mz];

                // Find the min/max indexes for the mz kernel.
                size_t i_min = 0;
                if (index_mz >= mz_kernel_hw) {
                    i_min = index_mz - mz_kernel_hw;
                }
                size_t i_max = grid.n - 1;
                if ((index_mz + mz_kernel_hw) < grid.n) {
                    i_max = index_mz + mz_kernel_hw;
                }

                for (size_t j = j_min; j <= j_max; ++j) {
                    for (size_t i = i_min; i <= i_max; ++i) {
                        double x = grid.bins_mz[i];
                        double y = grid.bins_rt[j];

                        // Calculate the Gaussian weight for this point.
                        double a = (x - current_mz) / sigma_mz;
                        double b = (y - current_rt) / sigma_rt;
                        double weight = std::exp(-0.5 * (a * a + b * b));

                        grid.data[i + j * n] += weight * current_intensity;
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
            grid.data[i] = grid.data[i] / weight;
        }
    }

    // Gaussian smoothing.
    //
    // The Gaussian 2D filter is separable. We obtain the same result with
    // faster performance by applying two 1D kernel convolutions instead. This
    // is specially noticeable on the full image.
    {
        auto smoothed_data = std::vector<double>(grid.n * grid.m);

        // Retention time smoothing.
        for (size_t j = 0; j < grid.m; ++j) {
            double current_rt = grid.bins_rt[j];
            size_t min_k = 0;
            if (j >= rt_kernel_hw) {
                min_k = j - rt_kernel_hw;
            }
            size_t max_k = grid.m - 1;
            if ((j + rt_kernel_hw) < grid.m) {
                max_k = j + rt_kernel_hw;
            }
            for (size_t i = 0; i < grid.n; ++i) {
                double sum_weights = 0;
                double sum_weighted_values = 0;
                for (size_t k = min_k; k <= max_k; ++k) {
                    double a = (current_rt - grid.bins_rt[k]) / sigma_rt;
                    double weight = std::exp(-0.5 * (a * a));
                    sum_weights += weight;
                    sum_weighted_values += weight * grid.data[i + k * grid.n];
                }
                smoothed_data[i + j * grid.n] =
                    sum_weighted_values / sum_weights;
            }
        }
        grid.data = smoothed_data;
    }
    {
        auto smoothed_data = std::vector<double>(grid.n * grid.m);

        // mz smoothing.
        //
        // Since sigma_mz is not constant, we need to calculate the kernel half
        // width for each potential mz value.
        for (size_t i = 0; i < grid.n; ++i) {
            double sigma_mz = sigma_mz_vec[i];
            double current_mz = grid.bins_mz[i];

            size_t min_k = 0;
            if (i >= mz_kernel_hw) {
                min_k = i - mz_kernel_hw;
            }
            size_t max_k = grid.n - 1;
            if ((i + mz_kernel_hw) < grid.n) {
                max_k = i + mz_kernel_hw;
            }
            for (size_t j = 0; j < grid.m; ++j) {
                double sum_weights = 0;
                double sum_weighted_values = 0;
                for (size_t k = min_k; k <= max_k; ++k) {
                    double a = (current_mz - grid.bins_mz[k]) / sigma_mz;
                    double weight = std::exp(-0.5 * (a * a));
                    sum_weights += weight;
                    sum_weighted_values += weight * grid.data[k + j * grid.n];
                }
                smoothed_data[i + j * grid.n] =
                    sum_weighted_values / sum_weights;
            }
        }
        grid.data = smoothed_data;
    }
    return grid;
}

Grid::Grid Grid::subset(Grid grid, double min_mz, double max_mz, double min_rt, double max_rt) {
    // Find min/max bin in mz and rt.
    size_t min_mz_idx = Search::lower_bound(grid.bins_mz, min_mz);
    size_t max_mz_idx = Search::lower_bound(grid.bins_mz, max_mz);
    size_t min_rt_idx = Search::lower_bound(grid.bins_rt, min_rt);
    size_t max_rt_idx = Search::lower_bound(grid.bins_rt, max_rt);

    // Initialize new grid.
    Grid new_grid;
    new_grid.n = max_mz_idx - min_mz_idx;
    new_grid.m = max_rt_idx - min_rt_idx;
    new_grid.k = grid.k;
    new_grid.t = grid.t;
    new_grid.instrument_type = grid.instrument_type;
    new_grid.reference_mz = grid.reference_mz;
    new_grid.fwhm_mz = grid.fwhm_mz;
    new_grid.fwhm_rt = grid.fwhm_rt;
    new_grid.min_mz = grid.bins_mz[min_mz_idx];
    new_grid.max_mz = grid.bins_mz[max_mz_idx];
    new_grid.min_rt = grid.bins_mz[min_rt_idx];
    new_grid.max_rt = grid.bins_mz[max_rt_idx];

    // Initialize new grid memory.
    new_grid.data = std::vector<double>(new_grid.n * new_grid.m);

    // Initialize bins.
    new_grid.bins_mz = std::vector<double>(new_grid.n);
    new_grid.bins_rt = std::vector<double>(new_grid.m);
    for (size_t i = 0; i < new_grid.n; i++) {
        size_t mz_idx = min_mz_idx + i;
        new_grid.bins_mz[i] = grid.bins_mz[mz_idx];
    }
    for (size_t j = 0; j < new_grid.m; j++) {
        size_t rt_idx = min_rt_idx + j;
        new_grid.bins_rt[j] = grid.bins_rt[rt_idx];
    }

    // Represent the bounds of this grid in mass-to-charge and retention time
    // coordinates.
    for (size_t j = 0; j < new_grid.m; j++) {
        for (size_t i = 0; i < new_grid.n; i++) {
            size_t mz_idx = min_mz_idx + i;
            size_t rt_idx = min_rt_idx + j;
            new_grid.data[i + j * new_grid.n] = grid.data[mz_idx + rt_idx * grid.n];
        }
    }

    return new_grid;
}
