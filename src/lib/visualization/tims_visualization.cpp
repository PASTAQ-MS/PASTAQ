#include "tims_visualization.h"
#include <cmath>

namespace timsvis
{

TimsVisualization::TimsVisualization(const timsdata::TimsData& tims_data)
    : tims_data_(tims_data) {}


HeatmapData TimsVisualization::heatmap_data(int x_bins, int y_bins) {
    double x_min = std::numeric_limits<double>::max();
    double x_max = std::numeric_limits<double>::lowest();
    double y_min = std::numeric_limits<double>::max();
    double y_max = std::numeric_limits<double>::lowest();

    // Retrieve all frames
    const auto& frames = tims_data_.getFrames();

    // Find min/max values
    for (const auto& frame_pair : frames) {
        const auto& frame = frame_pair.second;
        double y_value = frame.getTime(); // Assuming Time is the Y-axis

        if (y_value < y_min) y_min = y_value;
        if (y_value > y_max) y_max = y_value;

        for (const auto& spectrum : frame.getSpectra()) {
            for (double x_value : spectrum.getMz()) {
                if (x_value < x_min) x_min = x_value;
                if (x_value > x_max) x_max = x_value;
            }
        }
    }

    // Compute bin edges
    std::vector<double> x_edges(x_bins + 1);
    std::vector<double> y_edges(y_bins + 1);
    double x_bin_size = (x_max - x_min) / x_bins;
    double y_bin_size = (y_max - y_min) / y_bins;

    for (int i = 0; i <= x_bins; ++i) x_edges[i] = x_min + i * x_bin_size;
    for (int i = 0; i <= y_bins; ++i) y_edges[i] = y_min + i * y_bin_size;

    // Initialize matrix
    std::vector<std::vector<double>> value_matrix(x_bins, std::vector<double>(y_bins, 0.0));
    std::vector<std::vector<int>> count_matrix(x_bins, std::vector<int>(y_bins, 0));

    // Populate bins
    for (const auto& frame_pair : frames) {
        const auto& frame = frame_pair.second;
        int y_idx = std::min(static_cast<int>((frame.getTime() - y_min) / y_bin_size), y_bins - 1);

        for (const auto& spectrum : frame.getSpectra()) {
            const auto& x_values = spectrum.getMz();
            const auto& intensities = spectrum.getIntensity();
            size_t num_peaks = spectrum.getNumberOfPeaks();

            for (size_t i = 0; i < num_peaks; ++i) {
                double x_value = x_values[i];
                double value = intensities[i];

                int x_idx = std::min(static_cast<int>((x_value - x_min) / x_bin_size), x_bins - 1);

                value_matrix[x_idx][y_idx] += value;
                count_matrix[x_idx][y_idx] += 1;
            }
        }
    }

    // Compute log-transformed mean values
    for (int i = 0; i < x_bins; ++i) {
        for (int j = 0; j < y_bins; ++j) {
            if (count_matrix[i][j] > 0) {
                value_matrix[i][j] = std::log1p(value_matrix[i][j] / count_matrix[i][j]);
            }
        }
    }

    return {value_matrix, x_edges, y_edges};
}

}
