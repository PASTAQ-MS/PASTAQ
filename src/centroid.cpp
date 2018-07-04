#include "centroid.hpp"

std::vector<Centroid::Peak> Centroid::find_local_maxima(
    const Grid::Parameters &parameters, const std::vector<double> &data) {
    std::vector<Centroid::Peak> ret;
    unsigned int n_mz = parameters.dimensions.n;
    unsigned int n_rt = parameters.dimensions.m;
    for (size_t j = 1; j < n_rt; ++j) {
        for (size_t i = 1; i < n_mz; ++i) {
            int index = i + j * n_mz;

            double value = data[index];
            double right_value = data[index + 1];
            double left_value = data[index - 1];
            double top_value = data[index - n_mz];
            double top_right_value = data[index + 1 - n_mz];
            double top_left_value = data[index - 1 - n_mz];
            double bottom_value = data[index + n_mz];
            double bottom_right_value = data[index + 1 + n_mz];
            double bottom_left_value = data[index - 1 + n_mz];

            if ((value != 0) && (value > left_value) && (value > right_value) &&
                (value > top_value) && (value > top_right_value) &&
                (value > top_left_value) && (value > bottom_value) &&
                (value > bottom_left_value) && (value > bottom_right_value)) {
                Centroid::Peak peak = {};
                peak.i = i;
                peak.j = j;
                peak.height = value;
                ret.push_back(peak);
            }
        }
    }
    return ret;
}
