#include <thread>

#include "warp2d/warp2d_runners.hpp"

std::vector<Centroid::Peak> Warp2D::Runners::Parallel::run(
    const std::vector<Centroid::Peak>& target_peaks,
    const std::vector<Centroid::Peak>& source_peaks,
    const Warp2D::Parameters& parameters, uint64_t max_threads) {
    // Initialize parameters.
    int t = parameters.slack;        // Slack.
    int m = parameters.window_size;  // Segment/Window size.
    int nP = parameters.num_points;  // Number of points.
    int N = nP / m;                  // Number of segments.
    nP = N * m;

    // Find min/max retention times.
    double rt_min = std::numeric_limits<double>::infinity();
    double rt_max = -std::numeric_limits<double>::infinity();
    for (const auto& peak : target_peaks) {
        if (peak.rt < rt_min) {
            rt_min = peak.rt;
        }
        if (peak.rt > rt_max) {
            rt_max = peak.rt;
        }
    }
    for (const auto& peak : source_peaks) {
        if (peak.rt < rt_min) {
            rt_min = peak.rt;
        }
        if (peak.rt > rt_max) {
            rt_max = peak.rt;
        }
    }

    // TODO(alex): Is there a better way of doing this? I thought about adding
    // the equivalent of an extra sector at the beginning and end of the
    // retention time extremes, but maybe it's worse that this simple
    // implementation. Expand the time range by some factor, so we can have
    // warping at the peaks near the range limits (COW will not warp at the
    // range limits).
    double rt_expand_factor = 0.20;  // FIXME: Hardcoding this for now.
    rt_min -= (rt_max - rt_min) * rt_expand_factor;
    rt_max += (rt_max - rt_min) * rt_expand_factor;

    // The minimum time step.
    double delta_rt = (rt_max - rt_min) / (double)(nP - 1);
    double segment_rt_width = delta_rt * m;

    // Filter the peaks in each segment.
    int n_peaks_per_segment = 50;  // FIXME: Hardcoding this for now.
    std::vector<Centroid::Peak> target_peaks_filtered;
    std::vector<Centroid::Peak> source_peaks_filtered;
    for (int i = 0; i < N; ++i) {
        double rt_start = rt_min + i * segment_rt_width;
        double rt_end = rt_start + segment_rt_width;
        // Filter reference peaks.
        {
            auto peaks =
                Warp2D::peaks_in_rt_range(target_peaks, rt_start, rt_end);
            auto filtered_peaks = filter_peaks(peaks, n_peaks_per_segment);
            target_peaks_filtered.insert(end(target_peaks_filtered),
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

    // Prepare which group of levels we are going to send to every core. We
    // store the index of the levels into a groups array.
    auto groups = std::vector<std::vector<int>>(max_threads);
    size_t i = 0;
    size_t k = 0;
    while ((int)k < N) {
        groups[i].push_back(k);
        ++k;
        if (i == max_threads - 1) {
            i = 0;
        } else {
            ++i;
        }
    }

    std::vector<std::thread> threads(groups.size());
    for (size_t i = 0; i < groups.size(); ++i) {
        threads[i] = std::thread([i, &groups, &levels, rt_min, delta_rt,
                                  segment_rt_width, &target_peaks_filtered,
                                  &source_peaks_filtered]() {
            for (const auto& k : groups[i]) {
                auto& current_level = levels[k];
                double rt_start = rt_min + k * segment_rt_width;
                double rt_end = rt_start + segment_rt_width;
                Warp2D::compute_warped_similarities(
                    current_level, rt_start, rt_end, rt_min, delta_rt,
                    target_peaks_filtered, source_peaks_filtered);
            }
        });
    }
    for (auto& thread : threads) {
        thread.join();
    }

    auto warp_by = Warp2D::find_optimal_warping(levels, N);

    // Warp the sample peaks based on the optimal path.
    std::vector<Centroid::Peak> warped_peaks;
    warped_peaks.reserve(source_peaks.size());
    for (int i = 0; i < N; ++i) {
        double rt_start = rt_min + i * segment_rt_width;
        double rt_end = rt_start + segment_rt_width;

        int x_start = warp_by[i] + levels[i].start;
        int x_end = warp_by[i + 1] + levels[i + 1].start;

        double sample_rt_start = rt_min + x_start * delta_rt;
        double sample_rt_width = (x_end - x_start) * delta_rt;
        double sample_rt_end = sample_rt_start + sample_rt_width;

        auto warped_peaks_segment = Warp2D::warp_peaks(
            source_peaks, sample_rt_start, sample_rt_end, rt_start, rt_end);
        for (const auto& peak : warped_peaks_segment) {
            warped_peaks.push_back(peak);
        }
    }

    return warped_peaks;
}

std::vector<Centroid::Peak> Warp2D::Runners::Serial::run(
    const std::vector<Centroid::Peak>& target_peaks,
    const std::vector<Centroid::Peak>& source_peaks,
    const Warp2D::Parameters& parameters) {
    // Initialize parameters.
    int t = parameters.slack;        // Slack.
    int m = parameters.window_size;  // Segment/Window size.
    int nP = parameters.num_points;  // Number of points.
    int N = nP / m;                  // Number of segments.
    nP = N * m;

    // Find min/max retention times.
    double rt_min = std::numeric_limits<double>::infinity();
    double rt_max = -std::numeric_limits<double>::infinity();
    for (const auto& peak : target_peaks) {
        if (peak.rt < rt_min) {
            rt_min = peak.rt;
        }
        if (peak.rt > rt_max) {
            rt_max = peak.rt;
        }
    }
    for (const auto& peak : source_peaks) {
        if (peak.rt < rt_min) {
            rt_min = peak.rt;
        }
        if (peak.rt > rt_max) {
            rt_max = peak.rt;
        }
    }

    // TODO(alex): Is there a better way of doing this? I thought about adding
    // the equivalent of an extra sector at the beginning and end of the
    // retention time extremes, but maybe it's worse that this simple
    // implementation. Expand the time range by some factor, so we can have
    // warping at the peaks near the range limits (COW will not warp at the
    // range limits).
    double rt_expand_factor = 0.20;  // FIXME: Hardcoding this for now.
    rt_min -= (rt_max - rt_min) * rt_expand_factor;
    rt_max += (rt_max - rt_min) * rt_expand_factor;

    // The minimum time step.
    double delta_rt = (rt_max - rt_min) / (double)(nP - 1);
    double segment_rt_width = delta_rt * m;

    // Filter the peaks in each segment.
    int n_peaks_per_segment = 50;  // FIXME: Hardcoding this for now.
    std::vector<Centroid::Peak> target_peaks_filtered;
    std::vector<Centroid::Peak> source_peaks_filtered;
    for (int i = 0; i < N; ++i) {
        double rt_start = rt_min + i * segment_rt_width;
        double rt_end = rt_start + segment_rt_width;
        // Filter reference peaks.
        {
            auto peaks =
                Warp2D::peaks_in_rt_range(target_peaks, rt_start, rt_end);
            auto filtered_peaks = filter_peaks(peaks, n_peaks_per_segment);
            target_peaks_filtered.insert(end(target_peaks_filtered),
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

    for (int k = 0; k < N; ++k) {
        auto& current_level = levels[k];
        double rt_start = rt_min + k * segment_rt_width;
        double rt_end = rt_start + segment_rt_width;
        Warp2D::compute_warped_similarities(
            current_level, rt_start, rt_end, rt_min, delta_rt,
            target_peaks_filtered, source_peaks_filtered);
    }

    auto warp_by = Warp2D::find_optimal_warping(levels, N);

    // Warp the sample peaks based on the optimal path.
    std::vector<Centroid::Peak> warped_peaks;
    warped_peaks.reserve(source_peaks.size());
    for (int i = 0; i < N; ++i) {
        double rt_start = rt_min + i * segment_rt_width;
        double rt_end = rt_start + segment_rt_width;

        int x_start = warp_by[i] + levels[i].start;
        int x_end = warp_by[i + 1] + levels[i + 1].start;

        double sample_rt_start = rt_min + x_start * delta_rt;
        double sample_rt_width = (x_end - x_start) * delta_rt;
        double sample_rt_end = sample_rt_start + sample_rt_width;

        auto warped_peaks_segment = Warp2D::warp_peaks(
            source_peaks, sample_rt_start, sample_rt_end, rt_start, rt_end);
        for (const auto& peak : warped_peaks_segment) {
            warped_peaks.push_back(peak);
        }
    }

    return warped_peaks;
}
