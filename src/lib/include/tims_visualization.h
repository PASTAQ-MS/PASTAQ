#ifndef TIMS_VISUALIZATION_H
#define TIMS_VISUALIZATION_H

#include <vector>
#include "timsdatacpp.h"

using namespace timsdata;

namespace timsvis

{

struct HeatmapData {
    std::vector<std::vector<double>> value_matrix;
    std::vector<double> x_edges;
    std::vector<double> y_edges;
};

class TimsVisualization {
public:
    explicit TimsVisualization(const timsdata::TimsData& tims_data);
    HeatmapData heatmap_data(int x_bins, int y_bins);
private:
    const timsdata::TimsData& tims_data_; // Store reference
};

}

#endif // TIMS_VISUALIZATION_H
