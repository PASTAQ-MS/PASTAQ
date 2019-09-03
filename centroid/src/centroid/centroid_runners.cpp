#include <algorithm>
#include <thread>

#include "centroid_runners.hpp"

std::vector<Centroid::Peak> Centroid::Runners::Serial::run(
    const RawData::RawData &raw_data, const Grid::Mesh &mesh,
    size_t max_peaks) {
    // Finding local maxima.
    auto local_max = Centroid::find_local_maxima(mesh);

    // Sort the local_maxima by value.
    auto sort_local_max = [](const Centroid::LocalMax &p1,
                             const Centroid::LocalMax &p2) -> bool {
        return (p1.value >= p2.value);
    };
    std::stable_sort(local_max.begin(), local_max.end(), sort_local_max);

    std::vector<Centroid::Peak> peaks;
    for (const auto &max : local_max) {
        if (peaks.size() == max_peaks) {
            break;
        }
        auto peak = build_peak(raw_data, max);
        // FIXME: Number of raw points within the theoretical sigma
        // should be set by the user, with a sensible default. Same
        // with the minimum number of rt scans per peak.
        // Ensure peak quality.
        if (peak.raw_roi_num_points_within_sigma < 5 ||
            peak.raw_roi_num_scans < 3 || peak.raw_roi_sigma_mz == 0 ||
            peak.raw_roi_sigma_rt == 0) {
            continue;
        }
        peaks.push_back(peak);
    }

    // Update the peak ids.
    for (size_t i = 0; i < peaks.size(); ++i) {
        peaks[i].id = i;
    }

    return peaks;
}

std::vector<Centroid::Peak> Centroid::Runners::Parallel::run(
    const RawData::RawData &raw_data, const Grid::Mesh &mesh,
    size_t max_peaks) {
    // Finding local maxima.
    auto local_max = Centroid::find_local_maxima(mesh);

    // The number of groups/threads is set to the maximum possible concurrency.
    uint64_t max_threads = std::thread::hardware_concurrency();

    // Split the points into different groups for concurrency.
    std::vector<std::vector<size_t>> groups =
        std::vector<std::vector<size_t>>(max_threads);
    for (size_t i = 0; i < local_max.size(); ++i) {
        size_t k = i % max_threads;
        groups[k].push_back(i);
    }

    std::vector<std::thread> threads(max_threads);
    std::vector<std::vector<Centroid::Peak>> peaks_array(max_threads);
    for (size_t i = 0; i < groups.size(); ++i) {
        threads[i] = std::thread([&groups, &local_max, &peaks_array, &raw_data,
                                  &mesh, i]() {
            for (const auto &k : groups[i]) {
                auto peak = build_peak(raw_data, local_max[k]);
                // FIXME: Number of raw points within the theoretical sigma
                // should be set by the user, with a sensible default. Same
                // with the minimum number of rt scans per peak.
                // Ensure peak quality.
                if (peak.raw_roi_num_points_within_sigma < 5 ||
                    peak.raw_roi_num_scans < 3 || peak.raw_roi_sigma_mz == 0 ||
                    peak.raw_roi_sigma_rt == 0) {
                    continue;
                }
                peaks_array[i].push_back(peak);
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

    // Sort the peaks by height.
    auto sort_peaks = [](const Centroid::Peak &p1,
                         const Centroid::Peak &p2) -> bool {
        return (p1.local_max_height > p2.local_max_height) ||
               (p1.local_max_height == p2.local_max_height);
    };
    std::stable_sort(peaks.begin(), peaks.end(), sort_peaks);

    // Update the peak ids.
    for (size_t i = 0; i < peaks.size(); ++i) {
        peaks[i].id = i;
    }

    // Return maximum amount of peaks.
    // TODO: Figure a way of performing max peaks when multiple threads are
    // in place without having to go through all of them. Perhaps an atomic
    // operation for increment counter?
    if (peaks.size() > max_peaks) {
        peaks.resize(max_peaks);
    }

    return peaks;
}
