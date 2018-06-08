#include "tapp-mesh.hpp"

Mesh::Mesh(Grid::Dimensions dimensions, Grid::Bounds bounds)
    : m_data(dimensions.n * dimensions.m),
      m_dimensions(dimensions),
      m_bounds(bounds) {}

std::optional<double> Mesh::at(unsigned int i, unsigned int j) {
    if (m_data.size() == 0 || i > m_dimensions.n - 1 ||
        j > m_dimensions.m - 1) {
        return std::nullopt;
    }
    return m_data[i + j * m_dimensions.n];
}

bool Mesh::set(unsigned int i, unsigned int j, double value) {
    if (m_data.size() == 0 || i > m_dimensions.n - 1 ||
        j > m_dimensions.m - 1) {
        return false;
    }
    m_data[i + j * m_dimensions.n] = value;
    return true;
}

std::optional<double> Mesh::mz(unsigned int i) {
    if (m_data.size() == 0 || i > m_dimensions.n - 1) {
        return std::nullopt;
    }

    auto delta_mz = (m_bounds.max_mz - m_bounds.min_mz) /
                    static_cast<double>(m_dimensions.n - 1);
    return m_bounds.min_mz + delta_mz * i;
}

std::optional<double> Mesh::rt(unsigned int j) {
    if (m_data.size() == 0 || j > m_dimensions.m - 1) {
        return std::nullopt;
    }
    auto delta_rt = (m_bounds.max_rt - m_bounds.min_rt) /
                    static_cast<double>(m_dimensions.m - 1);
    return m_bounds.min_rt + delta_rt * j;
}

Grid::Dimensions Mesh::dim() { return m_dimensions; }

Grid::Bounds Mesh::bounds() { return m_bounds; }
