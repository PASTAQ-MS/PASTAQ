#include <cmath>
#include <iostream>

#include "tapp-weighted-mesh.hpp"

WeightedMesh::WeightedMesh(Grid::Dimensions dimensions, Grid::Bounds bounds)
    : Mesh(dimensions, bounds),
      m_weights(dimensions.n * dimensions.m),
      m_counts(dimensions.n * dimensions.m) {}

std::optional<double> WeightedMesh::weight_at(unsigned int i, unsigned int j) {
    if (m_data.size() == 0 || i > m_dimensions.n - 1 ||
        j > m_dimensions.m - 1) {
        return std::nullopt;
    }
    return m_weights[i + j * m_dimensions.m];
}

std::optional<double> WeightedMesh::counts_at(unsigned int i, unsigned int j) {
    if (m_data.size() == 0 || i > m_dimensions.n - 1 ||
        j > m_dimensions.m - 1) {
        return std::nullopt;
    }
    return m_counts[i + j * m_dimensions.m];
}

void WeightedMesh::print_all() {
    std::cout << "DATA:" << std::endl;
    for (unsigned int j = 0; j < m_dimensions.m; ++j) {
        for (unsigned int i = 0; i < m_dimensions.n; ++i) {
            auto e = m_data[i + j * m_dimensions.n];
            std::cout << e << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "WEIGHTS:" << std::endl;
    for (unsigned int j = 0; j < m_dimensions.m; ++j) {
        for (unsigned int i = 0; i < m_dimensions.n; ++i) {
            auto e = m_weights[i + j * m_dimensions.n];
            std::cout << e << "\t";
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
    std::cout << "COUNTS:" << std::endl;
    for (unsigned int j = 0; j < m_dimensions.m; ++j) {
        for (unsigned int i = 0; i < m_dimensions.n; ++i) {
            auto e = m_counts[i + j * m_dimensions.n];
            std::cout << e << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

std::optional<unsigned int> WeightedMesh::x_index(double mz) {
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

std::optional<unsigned int> WeightedMesh::y_index(double rt) {
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

bool WeightedMesh::splash(double mz_value, double rt_value, double value,
                          double sigma_mz, double sigma_rt) {
    // Get the gaussian square dimensions. Note that we are using a square
    // kernel approximation for simplicity and computational efficiency.
    double min_rt = rt_value - 2 * sigma_rt;
    double max_rt = rt_value + 2 * sigma_rt;
    double min_mz = mz_value - 2 * sigma_mz;
    double max_mz = mz_value + 2 * sigma_mz;

    // Even if the point lays outside the current grid, we still want to account
    // for it's contribution to the points in the frontier. However if the
    // minimum value in the splatting lays outside the grid, it will not have
    // any effect.
    auto i_min = x_index(min_mz);
    auto j_min = y_index(min_rt);
    if (i_min == std::nullopt && j_min == std::nullopt) {
        return false;
    }

    auto i_max = x_index(max_mz) ? x_index(max_mz) : m_dimensions.n - 1;
    auto j_max = y_index(max_rt) ? y_index(max_rt) : m_dimensions.m - 1;

    double x0 = mz_value;
    double y0 = rt_value;
    for (unsigned int j = j_min.value(); j <= j_max.value(); ++j) {
        for (unsigned int i = i_min.value(); i <= i_max.value(); ++i) {
            // No need to do boundary check, since we are sure we are inside the
            // grid.
            double x = mz(i).value();
            double y = rt(j).value();

            // Calculate the gaussian weight for this point.
            double a = (x - x0) / sigma_mz;
            double b = (y - y0) / sigma_rt;
            a *= a;
            b *= b;
            double weight = std::exp(-0.5 * (a + b));

            // Set the value, weight and counts.
            m_data[i + j * m_dimensions.n] += value * weight;
            m_weights[i + j * m_dimensions.n] += weight;
            ++m_counts[i + j * m_dimensions.n];
        }
    }
    return true;
}
