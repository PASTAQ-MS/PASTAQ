#ifndef TIMS_VISUALIZATION_H
#define TIMS_VISUALIZATION_H

#include <vector>
#include "timsdatacpp.h"

using namespace timsdata;

namespace timsvis

{

struct HeatmapData {
    std::vector<double> x_edges;
    std::vector<double> y_edges;
    std::vector<std::vector<double>> value_matrix;
};

// struct HeatmapData {
//     std::vector<double> x_values;  // X-axis bin edges
//     std::vector<double> y_values;  // Y-axis bin edges
//     std::vector<double> z_values;  // Intensity values (log-transformed)

//     // Optionally, you can add a constructor to make the initialization easier
//     HeatmapData(
//         const std::vector<double>& x_vals,
//         const std::vector<double>& y_vals,
//         const std::vector<double>& z_vals
//     ) : x_values(x_vals), y_values(y_vals), z_values(z_vals) {}
// };

class TimsVisualization {
public:
    explicit TimsVisualization(const timsdata::TimsData& tims_data);
    HeatmapData heatmap_data(int x_bins, int y_bins);
private:
    const timsdata::TimsData& tims_data_; // Store reference
};

}

#endif // TIMS_VISUALIZATION_H
