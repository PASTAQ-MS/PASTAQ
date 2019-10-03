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
uint64_t Grid::Mesh::y_index(double rt) {
    double delta_rt = this->fwhm_rt / this->k;
    return std::ceil((rt - this->min_rt) / delta_rt);
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
    uint64_t m = mesh.y_index(raw_data.max_rt) + 1;
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

    // Pre-calculate the kernel half widths for rt and mz.
    //
    // Since sigma_rt is constant, the size of the kernel will be the same
    // for the entire rt range. The same applies for mz, as we are using a
    // warped mesh that keeps the number of sampling points constant across the
    // entire m/z range.
    double delta_rt = mesh.fwhm_rt / num_samples_rt;
    double delta_mz = mesh.fwhm_mz / num_samples_mz;
    double sigma_mz_ref = RawData::fwhm_to_sigma(mesh.fwhm_mz);
    uint64_t rt_kernel_hw = 3 * sigma_rt / delta_rt;
    uint64_t mz_kernel_hw = 3 * sigma_mz_ref / delta_mz;

    // Gaussian splatting.
    {
        auto weights = std::vector<double>(n * m);
        for (size_t s = 0; s < raw_data.scans.size(); ++s) {
            const auto &scan = raw_data.scans[s];
            double current_rt = scan.retention_time;

            // Find the bin for the current retention time.
            size_t index_rt = mesh.y_index(current_rt);

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
                if (index_mz >= mz_kernel_hw) {
                    i_min = index_mz - mz_kernel_hw;
                }
                size_t i_max = mesh.n - 1;
                if ((index_mz + mz_kernel_hw) < mesh.n) {
                    i_max = index_mz + mz_kernel_hw;
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
            if (i >= mz_kernel_hw) {
                min_k = i - mz_kernel_hw;
            }
            size_t max_k = mesh.n - 1;
            if ((i + mz_kernel_hw) < mesh.n) {
                max_k = i + mz_kernel_hw;
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
