#include <cmath>

#include "tapp-mesh.hpp"

// TODO(alex): Handle case where bounds/dimensions are not correct.
RegularMesh::RegularMesh(Grid::Dimensions dimensions, Grid::Bounds bounds,
           Instrument::Type instrument_type,
           Grid::SmoothingParams smoothing_params)
    : m_data(dimensions.n * dimensions.m),
      m_dimensions(dimensions),
      m_bounds(bounds),
      m_instrument_type(instrument_type),
      m_smoothing_params(smoothing_params) {}

std::optional<double> RegularMesh::value_at(unsigned int i, unsigned int j) {
    if (m_data.empty() || i > m_dimensions.n - 1 || j > m_dimensions.m - 1) {
        return std::nullopt;
    }
    return m_data[i + j * m_dimensions.n];
}

bool RegularMesh::set_value(unsigned int i, unsigned int j, double value) {
    if (m_data.empty() || i > m_dimensions.n - 1 || j > m_dimensions.m - 1) {
        return false;
    }
    m_data[i + j * m_dimensions.n] = value;
    return true;
}

std::optional<double> RegularMesh::mz_at(unsigned int i) {
    if (m_data.empty() || i > m_dimensions.n - 1) {
        return std::nullopt;
    }

    auto delta_mz = (m_bounds.max_mz - m_bounds.min_mz) /
                    static_cast<double>(m_dimensions.n - 1);
    return m_bounds.min_mz + delta_mz * i;
}

std::optional<double> RegularMesh::rt_at(unsigned int j) {
    if (m_data.empty() || j > m_dimensions.m - 1) {
        return std::nullopt;
    }
    auto delta_rt = (m_bounds.max_rt - m_bounds.min_rt) /
                    static_cast<double>(m_dimensions.m - 1);
    return m_bounds.min_rt + delta_rt * j;
}

std::optional<unsigned int> RegularMesh::x_index(double mz) {
    // In order to be consistent, the maximum value is mz + delta_mz. This
    // ensures all intervals contain the same number of points.
    double delta_mz = (m_bounds.max_mz - m_bounds.min_mz) /
                      static_cast<double>(m_dimensions.n - 1);
    if (mz < m_bounds.min_mz || mz > m_bounds.max_mz + delta_mz) {
        return std::nullopt;
    }
    double d = mz - m_bounds.min_mz;
    auto i = static_cast<unsigned int>(d / delta_mz);
    if (i > m_dimensions.n - 1) {
        return std::nullopt;
    }
    return i;
}

std::optional<unsigned int> RegularMesh::y_index(double rt) {
    // In order to be consistent, the maximum value is rt + delta_rt. This
    // ensures all intervals contain the same number of points.
    double delta_rt = (m_bounds.max_rt - m_bounds.min_rt) /
                      static_cast<double>(m_dimensions.m - 1);
    if (rt < m_bounds.min_rt || rt > m_bounds.max_rt + delta_rt) {
        return std::nullopt;
    }
    double d = rt - m_bounds.min_rt;
    auto j = static_cast<unsigned int>(d / delta_rt);
    if (j > m_dimensions.m - 1) {
        return std::nullopt;
    }
    return j;
}

double RegularMesh::sigma_mz(double mz) {
    double sigma_mz = 0.0;
    switch (m_instrument_type) {
        case Instrument::ORBITRAP: {
            sigma_mz = m_smoothing_params.sigma_mz *
                       std::pow(mz / m_smoothing_params.mz, 1.5);
        } break;
        case Instrument::FTICR: {
            sigma_mz = m_smoothing_params.sigma_mz *
                       std::pow(mz / m_smoothing_params.mz, 2);
        } break;
        case Instrument::TOF: {
            sigma_mz = m_smoothing_params.sigma_mz * mz / m_smoothing_params.mz;
        } break;
        case Instrument::QUAD: {
            // QUAD/IONTRAP instruments maintain the same resolution across all
            // mass range.
            sigma_mz = m_smoothing_params.sigma_mz;
        } break;
    }
    return sigma_mz;
}

double RegularMesh::sigma_rt() { return m_smoothing_params.sigma_rt; }

Grid::Dimensions RegularMesh::dim() { return m_dimensions; }

Grid::Bounds RegularMesh::bounds() { return m_bounds; }
