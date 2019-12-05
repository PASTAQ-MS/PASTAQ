#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <thread>

#include "utils/interpolation.hpp"
#include "warp2d/warp2d.hpp"

std::vector<Centroid::Peak> Warp2D::peaks_in_rt_range(
    const std::vector<Centroid::Peak>& source_peaks, double time_start,
    double time_end) {
    std::vector<Centroid::Peak> ret;
    ret.reserve(source_peaks.size());
    size_t i = 0;
    for (const auto& peak : source_peaks) {
        if (peak.fitted_rt >= time_start && peak.fitted_rt < time_end) {
            ret.push_back(peak);
            ++i;
        }
    }
    ret.resize(i);
    return ret;
}

std::vector<Centroid::Peak> Warp2D::filter_peaks(
    std::vector<Centroid::Peak>& peaks, size_t n_peaks_max) {
    std::vector<Centroid::Peak> filtered_peaks;
    size_t n_peaks = n_peaks_max < peaks.size() ? n_peaks_max : peaks.size();
    if (n_peaks == 0) {
        return filtered_peaks;
    }
    filtered_peaks.reserve(n_peaks);
    auto sort_by_height = [](const Centroid::Peak& p1,
                             const Centroid::Peak& p2) -> bool {
        return (p1.fitted_height > p2.fitted_height);
    };
    std::sort(peaks.begin(), peaks.end(), sort_by_height);
    for (size_t i = 0; i < n_peaks; ++i) {
        auto peak = peaks[i];
        filtered_peaks.push_back(peak);
    }
    return filtered_peaks;
}

std::vector<Warp2D::Level> Warp2D::initialize_levels(int64_t N, int64_t m,
                                                     int64_t t, int64_t nP) {
    std::vector<Level> levels(N + 1);
    levels[N].start = nP;
    levels[N].end = nP;
    levels[N].nodes.push_back({0.0, 0});
    for (int64_t i = (N - 1); i >= 0; --i) {
        int64_t start = std::max((i * (m - t)), (nP - (N - i) * (m + t)));
        int64_t end = std::min((i * (m + t)), (nP - (N - i) * (m - t)));
        int64_t length = end - start + 1;
        levels[i].start = start;
        levels[i].end = end;
        levels[i].nodes = std::vector<Node>(length);
        for (int64_t j = 0; j < length; ++j) {
            levels[i].nodes[j].f = 0;
            levels[i].nodes[j].u =
                (levels[i + 1].end - levels[i + 1].start) / 2;
        }
    }

    // Calculate potential warpings on each level.
    for (int64_t k = 0; k < N; ++k) {
        auto& current_level = levels[k];
        const auto& next_level = levels[k + 1];

        for (int64_t i = 0; i < (int64_t)current_level.nodes.size(); ++i) {
            int64_t src_start = current_level.start + i;

            // The next node for the next level is subject to the following
            // constrains:
            //
            // x_{i + 1} = x_{i} + m + u, where u <- [-t, t]
            //
            int64_t x_end_min = std::max(src_start + m - t, next_level.start);
            int64_t x_end_max = std::min(src_start + m + t, next_level.end);
            int64_t j_min = x_end_min - next_level.start;
            int64_t j_max = x_end_max - next_level.start;
            for (int64_t j = j_min; j <= j_max; ++j) {
                int64_t src_end = next_level.start + j;
                levels[k].potential_warpings.push_back(
                    {i, j, src_start, src_end, 0});
            }
        }
    }
    return levels;
}

std::vector<Centroid::Peak> Warp2D::interpolate_peaks(
    const std::vector<Centroid::Peak>& source_peaks, double source_rt_start,
    double source_rt_end, double ref_rt_start, double ref_rt_end) {
    auto warped_peaks =
        Warp2D::peaks_in_rt_range(source_peaks, source_rt_start, source_rt_end);

    for (auto& peak : warped_peaks) {
        double previous_rt = peak.fitted_rt;
        double x =
            (previous_rt - source_rt_start) / (source_rt_end - source_rt_start);
        double new_rt = Interpolation::lerp(ref_rt_start, ref_rt_end, x);
        peak.rt_delta = new_rt - previous_rt;
    }
    return warped_peaks;
}

void Warp2D::compute_warped_similarities(
    Warp2D::Level& level, double rt_start, double rt_end, double rt_min,
    double delta_rt, const std::vector<Centroid::Peak>& ref_peaks,
    const std::vector<Centroid::Peak>& source_peaks) {
    auto ref_peaks_segment =
        Warp2D::peaks_in_rt_range(ref_peaks, rt_start, rt_end);

    for (auto& warping : level.potential_warpings) {
        int64_t src_start = warping.src_start;
        int64_t src_end = warping.src_end;

        double sample_rt_start = rt_min + warping.src_start * delta_rt;
        double sample_rt_width = (src_end - src_start) * delta_rt;
        double sample_rt_end = sample_rt_start + sample_rt_width;

        auto source_peaks_warped = Warp2D::interpolate_peaks(
            source_peaks, sample_rt_start, sample_rt_end, rt_start, rt_end);

        double similarity = Centroid::cumulative_overlap(ref_peaks_segment,
                                                         source_peaks_warped);
        warping.warped_similarity = similarity;
    }
}

std::vector<int64_t> Warp2D::find_optimal_warping(std::vector<Level>& levels) {
    int64_t N = levels.size() - 1;

    // Perform dynamic programming to update the cumulative similarity and the
    // optimal predecessor.
    for (int64_t k = N - 1; k >= 0; --k) {
        auto& current_level = levels[k];
        const auto& next_level = levels[k + 1];
        // TODO: If all nodes within the posibilities for this anchor have
        // the same performance (Below a given threshold), DO NOT MOVE THE
        // ANCHOR.
        for (const auto& warping : current_level.potential_warpings) {
            auto& node_i = current_level.nodes[warping.i];
            const auto& node_j = next_level.nodes[warping.j];
            double f_sum = node_j.f + warping.warped_similarity;
            if (f_sum > node_i.f) {
                node_i.f = f_sum;
                node_i.u = warping.j;
            }
        }
    }

    // Walk forward nodes to find optimal warping path.
    std::vector<int64_t> warp_by;
    warp_by.reserve(N + 1);
    warp_by.push_back(0);
    for (int64_t i = 0; i < N; ++i) {
        auto u = levels[i].nodes[warp_by[i]].u;
        warp_by.push_back(u);
    }
    return warp_by;
}

std::vector<Centroid::Peak> Warp2D::warp_peaks(
    const std::vector<Centroid::Peak>& source_peaks,
    const Warp2D::TimeMap& time_map) {
    std::vector<Centroid::Peak> warped_peaks;
    warped_peaks.reserve(source_peaks.size());
    for (size_t i = 0; i < time_map.num_segments; ++i) {
        double rt_start = time_map.rt_start[i];
        double rt_end = time_map.rt_end[i];
        double sample_rt_start = time_map.sample_rt_start[i];
        double sample_rt_end = time_map.sample_rt_end[i];
        auto warped_peaks_segment = Warp2D::interpolate_peaks(
            source_peaks, sample_rt_start, sample_rt_end, rt_start, rt_end);
        for (const auto& peak : warped_peaks_segment) {
            warped_peaks.push_back(peak);
        }
    }
    std::sort(warped_peaks.begin(), warped_peaks.end(),
              [](const Centroid::Peak& p1, const Centroid::Peak& p2) -> bool {
                  return (p1.id < p2.id);
              });
    return warped_peaks;
}

Warp2D::TimeMap Warp2D::calculate_time_map(
    const std::vector<Centroid::Peak>& ref_peaks,
    const std::vector<Centroid::Peak>& source_peaks,
    const Warp2D::Parameters& parameters, uint64_t max_threads) {
    // Initialize parameters.
    int n_peaks_per_segment =
        parameters.peaks_per_window;  // Maximum number of peaks on a window.
    int t = parameters.slack;         // Slack.
    int m = parameters.window_size;   // Segment/Window size.
    int nP = parameters.num_points;   // Number of points.
    int N = nP / m;                   // Number of segments.
    nP = N * m;

    // Find min/max retention times.
    double rt_min = std::numeric_limits<double>::infinity();
    double rt_max = -std::numeric_limits<double>::infinity();
    for (const auto& peak : ref_peaks) {
        if (peak.fitted_rt < rt_min) {
            rt_min = peak.fitted_rt;
        }
        if (peak.fitted_rt > rt_max) {
            rt_max = peak.fitted_rt;
        }
    }
    for (const auto& peak : source_peaks) {
        if (peak.fitted_rt < rt_min) {
            rt_min = peak.fitted_rt;
        }
        if (peak.fitted_rt > rt_max) {
            rt_max = peak.fitted_rt;
        }
    }

    rt_min -= (rt_max - rt_min) * parameters.rt_expand_factor;
    rt_max += (rt_max - rt_min) * parameters.rt_expand_factor;

    // The minimum time step.
    double delta_rt = (rt_max - rt_min) / (double)(nP - 1);
    double segment_rt_width = delta_rt * m;

    // Filter the peaks in each segment.
    std::vector<Centroid::Peak> ref_peaks_filtered;
    std::vector<Centroid::Peak> source_peaks_filtered;
    for (int i = 0; i < N; ++i) {
        double rt_start = rt_min + i * segment_rt_width;
        double rt_end = rt_start + segment_rt_width;
        // Filter reference peaks.
        {
            auto peaks = Warp2D::peaks_in_rt_range(ref_peaks, rt_start, rt_end);
            auto filtered_peaks = filter_peaks(peaks, n_peaks_per_segment);
            ref_peaks_filtered.insert(end(ref_peaks_filtered),
                                      begin(filtered_peaks),
                                      end(filtered_peaks));
        }
        // Filter source peaks.
        {
            auto peaks =
                Warp2D::peaks_in_rt_range(source_peaks, rt_start, rt_end);
            auto filtered_peaks = filter_peaks(peaks, n_peaks_per_segment);
            source_peaks_filtered.insert(end(source_peaks_filtered),
                                         begin(filtered_peaks),
                                         end(filtered_peaks));
        }
    }

    // Initialize nodes.
    auto levels = Warp2D::initialize_levels(N, m, t, nP);

    // Prepare maximum concurrency.
    uint64_t num_threads = std::thread::hardware_concurrency();
    if (num_threads > max_threads) {
        num_threads = max_threads;
    }

    // Prepare which group of levels we are going to send to every core. We
    // store the index of the levels into a groups array.
    auto groups = std::vector<std::vector<int>>(num_threads);
    size_t i = 0;
    size_t k = 0;
    while ((int)k < N) {
        groups[i].push_back(k);
        ++k;
        if (i == num_threads - 1) {
            i = 0;
        } else {
            ++i;
        }
    }

    std::vector<std::thread> threads(groups.size());
    for (size_t i = 0; i < groups.size(); ++i) {
        threads[i] = std::thread([i, &groups, &levels, rt_min, delta_rt,
                                  segment_rt_width, &ref_peaks_filtered,
                                  &source_peaks_filtered]() {
            for (const auto& k : groups[i]) {
                auto& current_level = levels[k];
                double rt_start = rt_min + k * segment_rt_width;
                double rt_end = rt_start + segment_rt_width;
                Warp2D::compute_warped_similarities(
                    current_level, rt_start, rt_end, rt_min, delta_rt,
                    ref_peaks_filtered, source_peaks_filtered);
            }
        });
    }
    for (auto& thread : threads) {
        thread.join();
    }

    auto warp_by = Warp2D::find_optimal_warping(levels);

    // Build the TimeMap.
    TimeMap time_map = {};
    time_map.rt_min = rt_min;
    time_map.rt_max = rt_max;
    time_map.num_segments = N;
    for (int i = 0; i < N; ++i) {
        double rt_start = rt_min + i * segment_rt_width;
        double rt_end = rt_start + segment_rt_width;

        int x_start = warp_by[i] + levels[i].start;
        int x_end = warp_by[i + 1] + levels[i + 1].start;

        double sample_rt_start = rt_min + x_start * delta_rt;
        double sample_rt_width = (x_end - x_start) * delta_rt;
        double sample_rt_end = sample_rt_start + sample_rt_width;

        time_map.rt_start.push_back(rt_start);
        time_map.rt_end.push_back(rt_end);
        time_map.sample_rt_start.push_back(sample_rt_start);
        time_map.sample_rt_end.push_back(sample_rt_end);
    }

    return time_map;
}

double Warp2D::warp(const Warp2D::TimeMap& time_map, double rt) {
    // Find the segment that contains the given rt.
    uint64_t segment = 0;
    for (size_t i = 0; i < time_map.num_segments; ++i) {
        if (rt >= time_map.sample_rt_start[i] &&
            rt < time_map.sample_rt_end[i]) {
            segment = i;
            break;
        }
    }
    // Interpolate.
    double rt_start = time_map.rt_start[segment];  // After warping
    double rt_end = time_map.rt_end[segment];      // After warping
    double sample_rt_start = time_map.sample_rt_start[segment];  // Original
    double sample_rt_end = time_map.sample_rt_end[segment];      // Original
    double x = (rt - sample_rt_start) / (sample_rt_end - sample_rt_start);
    return Interpolation::lerp(rt_start, rt_end, x);
}
