#include <algorithm>
#include <cmath>

#include "centroid/centroid.hpp"

std::vector<Centroid::Point> Centroid::find_local_maxima(
    const Centroid::Parameters &parameters, const std::vector<double> &data) {
    std::vector<Centroid::Point> points;
    auto grid_params = parameters.grid_params;
    uint64_t n_mz = grid_params.dimensions.n;
    uint64_t n_rt = grid_params.dimensions.m;
    for (size_t j = 1; j < n_rt; ++j) {
        for (size_t i = 1; i < n_mz; ++i) {
            int index = i + j * n_mz;

            // NOTE(alex): The definition of a local maxima in a 2D space might
            // have different interpretations. i.e. We can select the 8
            // neighbours and the local maxima will be marked if all points are
            // below the central value. Alternatively, only a number N of
            // neighbours can be used, for example only the 4 cardinal
            // directions from the value under study.
            //
            // ----------------------------------------------
            // |              | top_value    |              |
            // ----------------------------------------------
            // | left_value   | value        | right_value  |
            // ----------------------------------------------
            // |              | bottom_value |              |
            // ----------------------------------------------
            double value = data[index];
            double right_value = data[index + 1];
            double left_value = data[index - 1];
            double top_value = data[index - n_mz];
            double bottom_value = data[index + n_mz];

            if ((value != 0) && (value > left_value) && (value > right_value) &&
                (value > top_value) && (value > bottom_value)) {
                Centroid::Point point = {};
                point.i = i;
                point.j = j;
                point.value = value;
                points.push_back(point);
            }
        }
    }

    // Sort the local maxima by descending height.
    auto sort_by_value = [](const Centroid::Point &p1,
                            const Centroid::Point &p2) -> bool {
        return (p1.value > p2.value);
    };
    std::stable_sort(points.begin(), points.end(), sort_by_value);

    if (parameters.n_peaks != 0 && parameters.n_peaks < points.size()) {
        points.resize(parameters.n_peaks);
    }
    return points;
}

// FIXME: Remove this function in favour of the new build_peak function.
std::optional<Centroid::Peak> Centroid::build_peak(
    const Centroid::Point &local_max, const Centroid::Parameters &parameters,
    const std::vector<double> &data) {
    Centroid::Peak peak = {};
    return peak;
}

std::tuple<std::vector<double>, std::vector<double>> Centroid::Peak::xic(
    const RawData::RawData &raw_data, std::string method) {
    std::vector<double> rt;
    std::vector<double> intensity;
    Centroid::Peak &peak = *this;
    return raw_data.xic(peak.roi_min_mz, peak.roi_max_mz, peak.roi_min_rt,
                        peak.roi_max_rt, method);
}
