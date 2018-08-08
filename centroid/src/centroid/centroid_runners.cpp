#include <thread>

#include "centroid_runners.hpp"

std::vector<Centroid::Peak> Centroid::Runners::Serial::run(
    const Centroid::Parameters &parameters, const std::vector<double> &data) {
    // Finding local maxima.
    auto local_max_points = Centroid::find_local_maxima(parameters, data);

    // Building peaks.
    std::vector<Centroid::Peak> peaks;
    for (const auto &point : local_max_points) {
        auto peak = Centroid::build_peak(point, parameters, data);
        if (peak) {
            peaks.push_back(peak.value());
        }
    }
    return peaks;
}

std::vector<Centroid::Peak> Centroid::Runners::Parallel::run(
    uint64_t max_threads, const Centroid::Parameters &parameters,
    const std::vector<double> &data) {
    // Finding local maxima.
    auto local_max_points = Centroid::find_local_maxima(parameters, data);

    // Split the points into different groups for concurrency.
    uint64_t n_points;
    if (local_max_points.size() % max_threads == 0) {
        n_points = local_max_points.size() / max_threads;
    } else {
        n_points = local_max_points.size() / (max_threads - 1);
    }
    std::vector<std::vector<Centroid::Point>> groups;
    for (size_t i = 0; i < max_threads; ++i) {
        auto begin = i * n_points;
        auto end = i * n_points + n_points;
        if (end > local_max_points.size()) {
            end = local_max_points.size();
        }
        if (begin == end) {
            break;
        }
        std::vector<Centroid::Point> group_points(&local_max_points[begin],
                                                  &local_max_points[end]);
        groups.push_back(group_points);
    }

    std::vector<std::thread> threads(groups.size());
    std::vector<std::vector<Centroid::Peak>> peaks_array(groups.size());
    for (size_t i = 0; i < groups.size(); ++i) {
        threads[i] =
            std::thread([&groups, &parameters, &peaks_array, &data, i]() {
                // Building peaks.
                for (const auto &point : groups[i]) {
                    auto peak = Centroid::build_peak(point, parameters, data);
                    if (peak) {
                        peaks_array[i].push_back(peak.value());
                    }
                }
            });
    }

    // Wait for the threads to finish.
    for (auto &thread : threads) {
        thread.join();
    }

    // Join peak groups.
    std::vector<Centroid::Peak> peaks;
    for (size_t i = 0; i < peaks_array.size(); ++i) {
        peaks.insert(end(peaks), begin(peaks_array[i]), end(peaks_array[i]));
    }
    return peaks;
}
