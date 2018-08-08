#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <thread>

#include "warp2d/warp2d.hpp"

// TODO(alex): move to header?
struct Node {
    double f;  // Cummulative similarity value.
    int u;     // Optimal predecessor index.
};

struct PotentialWarping {
    int i;  // Index of the node on the current level for x_start.
    int j;  // Index of the node on the next level for x_end.
    int x_start;
    int x_end;
    double warped_similarity;
};

// TODO(alex): move to header?
struct Level {
    int start;
    int end;
    std::vector<Node> nodes;
    std::vector<PotentialWarping> potential_warpings;
};

std::vector<Centroid::Peak> peaks_in_rt_range(
    const std::vector<Centroid::Peak>& source_peaks, double time_start,
    double time_end) {
    std::vector<Centroid::Peak> ret;
    ret.reserve(source_peaks.size());
    int i = 0;
    for (const auto& peak : source_peaks) {
        if (peak.rt >= time_start && peak.rt < time_end) {
            ret.push_back(peak);
            ++i;
        }
        // TODO(alex): Should we use rt or rt_centroid?
    }
    ret.resize(i);
    return ret;
}

// TODO(alex): Move to utils module.
// Perform numerically stable linear interpolation of x between y_0 and y_1.
// Note that x is a number between 0 and 1. This is the equivalent of the
// following formula:
//
//     (y - y_0) / (y_1 - y_0) = x;
//     y = x * (y_1 - y_0) + y_0;
//
double lerp(double y_0, double y_1, double x) {
    return (1 - x) * y_0 + x * y_1;
}

double Warp2D::peak_overlap(const Centroid::Peak& peak_a,
                            const Centroid::Peak& peak_b) {
    // TODO(alex): Use rt/mz/height or rt_centroid/mz_centroid/height_centroid?

    // Early return if the peaks do not intersect in the +/-2 * sigma_mz/rt
    {
        double min_rt_a = peak_a.rt - 3 * peak_a.sigma_rt;
        double max_rt_a = peak_a.rt + 3 * peak_a.sigma_rt;
        double min_mz_a = peak_a.mz - 3 * peak_a.sigma_mz;
        double max_mz_a = peak_a.mz + 3 * peak_a.sigma_mz;
        double min_rt_b = peak_b.rt - 3 * peak_b.sigma_rt;
        double max_rt_b = peak_b.rt + 3 * peak_b.sigma_rt;
        double min_mz_b = peak_b.mz - 3 * peak_b.sigma_mz;
        double max_mz_b = peak_b.mz + 3 * peak_b.sigma_mz;

        if (max_rt_a < min_rt_b || max_rt_b < min_rt_a || max_mz_a < min_mz_b ||
            max_mz_b < min_mz_a) {
            return 0;
        }
    }

    // Calculate the gaussian contribution of the overlap between two points in
    // one dimension.
    auto gaussian_contribution = [](double x_a, double x_b, double sigma_a,
                                    double sigma_b) -> double {
        double var_a = std::pow(sigma_a, 2);
        double var_b = std::pow(sigma_b, 2);

        double a = (var_a + var_b) / (var_a * var_b) *
                   std::pow((x_a * var_b + x_b * var_a) / (var_a + var_b), 2);
        double b = (x_a * x_a) / var_a + (x_b * x_b) / var_b;

        return std::exp(0.5 * (a - b)) / std::sqrt(var_a + var_b);
    };

    auto rt_contrib = gaussian_contribution(peak_a.rt, peak_b.rt,
                                            peak_a.sigma_rt, peak_b.sigma_rt);
    auto mz_contrib = gaussian_contribution(peak_a.mz, peak_b.mz,
                                            peak_a.sigma_mz, peak_b.sigma_mz);

    return rt_contrib * mz_contrib * peak_a.height * peak_b.height;
}

double Warp2D::similarity_2D(const std::vector<Centroid::Peak>& set_a,
                             const std::vector<Centroid::Peak>& set_b) {
    double cummulative_similarity = 0;
    for (const auto& peak_a : set_a) {
        for (const auto& peak_b : set_b) {
            cummulative_similarity += Warp2D::peak_overlap(peak_a, peak_b);
        }
    }
    return cummulative_similarity;
}

// Filter the peaks based on peak height. Note that this function modifies the
// given `peaks` argument by sorting the vector in place.
std::vector<Centroid::Peak> filter_peaks(std::vector<Centroid::Peak>& peaks,
                                         size_t n_peaks_max) {
    std::vector<Centroid::Peak> filtered_peaks;
    size_t n_peaks = n_peaks_max < peaks.size() ? n_peaks_max : peaks.size();
    if (n_peaks == 0) {
        return filtered_peaks;
    }
    filtered_peaks.reserve(n_peaks);
    auto sort_by_height = [](const Centroid::Peak& p1,
                             const Centroid::Peak& p2) -> bool {
        return (p1.height > p2.height);
    };
    std::sort(peaks.begin(), peaks.end(), sort_by_height);
    for (size_t i = 0; i < n_peaks; ++i) {
        auto peak = peaks[i];
        // We don't need the inner points or boundary points for this algorithm.
        // In order to reduce the number of allocations and memory consumption
        // we remove these from the peak object.
        peak.points = {};
        peak.boundary = {};
        filtered_peaks.push_back(peak);
    }
    return filtered_peaks;
}

// TODO(alex): What do we do with rt_centroid? Currently we are only warping
// Centroid::Peak.rt...
// NOTE: Serial version.
// std::vector<Centroid::Peak> Warp2D::warp_peaks(
// const std::vector<Centroid::Peak>& target_peaks,
// const std::vector<Centroid::Peak>& source_peaks,
// const Warp2D::Parameters& parameters) {
//// Initialize parameters.
// int t = parameters.slack;        // Slack.
// int m = parameters.window_size;  // Segment/Window size.
// int nP = parameters.num_points;  // Number of points.
// int N = nP / m;                  // Number of segments.
// nP = N * m;

//// Find min/max retention times.
// double rt_min = std::numeric_limits<double>::infinity();
// double rt_max = -std::numeric_limits<double>::infinity();
// for (const auto& peak : target_peaks) {
// if (peak.rt < rt_min) {
// rt_min = peak.rt;
//}
// if (peak.rt > rt_max) {
// rt_max = peak.rt;
//}
//}
// for (const auto& peak : source_peaks) {
// if (peak.rt < rt_min) {
// rt_min = peak.rt;
//}
// if (peak.rt > rt_max) {
// rt_max = peak.rt;
//}
//}

//// TODO(alex): Is there a better way of doing this? I thought about adding
//// the equivalent of an extra sector at the beginning and end of the
//// retention time extremes, but maybe it's worse that this simple
//// implementation. Expand the time range by some factor, so we can have
//// warping at the peaks near the range limits (COW will not warp at the
//// range limits).
// double rt_expand_factor = 0.20;  // FIXME: Hardcoding this for now.
// rt_min -= (rt_max - rt_min) * rt_expand_factor;
// rt_max += (rt_max - rt_min) * rt_expand_factor;

//// The minimum time step.
// double delta_rt = (rt_max - rt_min) / (double)(nP - 1);
// double segment_rt_width = delta_rt * m;

//// Filter the peaks in each segment.
// int n_peaks_per_segment = 50;  // FIXME: Hardcoding this for now.
// std::vector<Centroid::Peak> target_peaks_filtered;
// std::vector<Centroid::Peak> source_peaks_filtered;
// for (int i = 0; i < N; ++i) {
// double rt_start = rt_min + i * segment_rt_width;
// double rt_end = rt_start + segment_rt_width;
//// Filter reference peaks.
//{
// auto peaks = peaks_in_rt_range(target_peaks, rt_start, rt_end);
// auto filtered_peaks = filter_peaks(peaks, n_peaks_per_segment);
// target_peaks_filtered.insert(end(target_peaks_filtered),
// begin(filtered_peaks),
// end(filtered_peaks));
//}
//// Filter source peaks.
//{
// auto peaks = peaks_in_rt_range(source_peaks, rt_start, rt_end);
// auto filtered_peaks = filter_peaks(peaks, n_peaks_per_segment);
// source_peaks_filtered.insert(end(source_peaks_filtered),
// begin(filtered_peaks),
// end(filtered_peaks));
//}
//}

//// Initialize nodes.
// std::vector<Level> levels(N + 1);
// for (int i = 0; i < N; ++i) {
// int start = std::max((i * (m - t)), (nP - (N - i) * (m + t)));
// int end = std::min((i * (m + t)), (nP - (N - i) * (m - t)));
// int length = end - start + 1;
// levels[i].start = start;
// levels[i].end = end;
// levels[i].nodes = std::vector<Node>(length);
// for (int j = 0; j < length; ++j) {
// levels[i].nodes[j].f = -std::numeric_limits<double>::infinity();
// levels[i].nodes[j].u = 0;
//}
//}
// levels[N].start = nP;
// levels[N].end = nP;
// levels[N].nodes.push_back({0.0, 0});

//// Perform dynamic programming to find optimal warping path.
// for (int i = N - 1; i >= 0; --i) {
// auto& current_level = levels[i];
// const auto& next_level = levels[i + 1];

// double rt_start = rt_min + i * segment_rt_width;
// double rt_end = rt_start + segment_rt_width;
// auto target_peaks_segment =
// peaks_in_rt_range(target_peaks_filtered, rt_start, rt_end);

//// We want to update the nodes at the current level based on the optimal
//// combination of warping between the current (node_i) and next level
//// nodes (node_j).
// for (size_t k = 0; k < current_level.nodes.size(); ++k) {
// int x_start = current_level.start + k;
// auto& node_i = current_level.nodes[k];

//// The next node for the next level is subject to the following
//// constrains:
////
//// x_{i + 1} = x_{i} + m + u, where u <- [-t, t]
////
// int x_end_min = std::max(x_start + m - t, next_level.start);
// int x_end_max = std::min(x_start + m + t, next_level.end);
// int j_min = x_end_min - next_level.start;
// int j_max = x_end_max - next_level.start;
// for (int j = j_min; j <= j_max; ++j) {
// int x_end = next_level.start + j;
// auto& node_j = next_level.nodes[j];

// double sample_rt_start = rt_min + x_start * delta_rt;
// double sample_rt_width = (x_end - x_start) * delta_rt;
// double sample_rt_end = sample_rt_start + sample_rt_width;

// auto source_peaks_warped = peaks_in_rt_range(
// source_peaks_filtered, sample_rt_start, sample_rt_end);

//// Warp the peaks by linearly interpolating their retention time
//// to the current segment's. Note that we are just performing
//// linear displacement of the center of the peaks, we do not
//// deform the peak shape by adjusting the sigmas.
// for (auto& peak : source_peaks_warped) {
// double x = (peak.rt - sample_rt_start) / sample_rt_width;
// peak.rt = lerp(rt_start, rt_end, x);
//}

//// Calculate the peak overlap between the reference and warped
//// peaks for this segment.
// double similarity = Warp2D::similarity_2D(target_peaks_segment,
// source_peaks_warped);
// double f_sum = node_j.f + similarity;
// if (f_sum > node_i.f) {
// node_i.f = f_sum;
// node_i.u = j;
//}
//}
//}
//}

//// Walk back nodes to find optimal warping path.
// std::vector<int> warp_by;
// warp_by.reserve(N + 1);
// warp_by.push_back(0);
// for (int i = 0; i < N; ++i) {
// auto u = levels[i].nodes[warp_by[i]].u;
// warp_by.push_back(u);
//}

//// Warp the sample peaks based on the optimal path.
// std::vector<Centroid::Peak> warped_peaks;
// warped_peaks.reserve(source_peaks.size());
// for (int i = 0; i < N; ++i) {
// double rt_start = rt_min + i * segment_rt_width;
// double rt_end = rt_start + segment_rt_width;

// int x_start = warp_by[i] + levels[i].start;
// int x_end = warp_by[i + 1] + levels[i + 1].start;

// double sample_rt_start = rt_min + x_start * delta_rt;
// double sample_rt_width = (x_end - x_start) * delta_rt;
// double sample_rt_end = sample_rt_start + sample_rt_width;
// auto source_peaks_segment =
// peaks_in_rt_range(source_peaks, sample_rt_start, sample_rt_end);

//// Make a copy of the peaks for warping.
// std::vector<Centroid::Peak> source_peaks_warped;
// source_peaks_warped.reserve(source_peaks_segment.size());
// for (const auto& peak : source_peaks_segment) {
// source_peaks_warped.push_back(peak);
//}

//// Warp the peaks.
// for (auto& peak : source_peaks_warped) {
// double x = (peak.rt - sample_rt_start) / sample_rt_width;
// peak.rt = lerp(rt_start, rt_end, x);
// warped_peaks.push_back(peak);
//}
//}

// return warped_peaks;
//}

// NOTE: Parallel version
std::vector<Centroid::Peak> Warp2D::warp_peaks(
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
            auto peaks = peaks_in_rt_range(target_peaks, rt_start, rt_end);
            auto filtered_peaks = filter_peaks(peaks, n_peaks_per_segment);
            target_peaks_filtered.insert(end(target_peaks_filtered),
                                         begin(filtered_peaks),
                                         end(filtered_peaks));
        }
        // Filter source peaks.
        {
            auto peaks = peaks_in_rt_range(source_peaks, rt_start, rt_end);
            auto filtered_peaks = filter_peaks(peaks, n_peaks_per_segment);
            source_peaks_filtered.insert(end(source_peaks_filtered),
                                         begin(filtered_peaks),
                                         end(filtered_peaks));
        }
    }

    // Initialize nodes.
    std::vector<Level> levels(N + 1);
    for (int i = 0; i < N; ++i) {
        int start = std::max((i * (m - t)), (nP - (N - i) * (m + t)));
        int end = std::min((i * (m + t)), (nP - (N - i) * (m - t)));
        int length = end - start + 1;
        levels[i].start = start;
        levels[i].end = end;
        levels[i].nodes = std::vector<Node>(length);
        for (int j = 0; j < length; ++j) {
            levels[i].nodes[j].f = -std::numeric_limits<double>::infinity();
            levels[i].nodes[j].u = 0;
        }
    }
    levels[N].start = nP;
    levels[N].end = nP;
    levels[N].nodes.push_back({0.0, 0});

    // Calculate potential warpings on each level.
    for (int k = 0; k < N; ++k) {
        auto& current_level = levels[k];
        const auto& next_level = levels[k + 1];

        for (int i = 0; i < (int)current_level.nodes.size(); ++i) {
            int x_start = current_level.start + i;

            // The next node for the next level is subject to the following
            // constrains:
            //
            // x_{i + 1} = x_{i} + m + u, where u <- [-t, t]
            //
            int x_end_min = std::max(x_start + m - t, next_level.start);
            int x_end_max = std::min(x_start + m + t, next_level.end);
            int j_min = x_end_min - next_level.start;
            int j_max = x_end_max - next_level.start;
            for (int j = j_min; j <= j_max; ++j) {
                int x_end = next_level.start + j;
                levels[k].potential_warpings.push_back(
                    {i, j, x_start, x_end, 0});
            }
        }
    }

    // Prepare which group of levels we are going to send to every core. We
    // store the index of the levels into a groups array.
    std::vector<std::vector<int>> groups = {};
    int max_threads =
        std::thread::hardware_concurrency();  // FIXME: Hardcoding this for now
    std::cout << "distributing groups..." << std::endl;
    if (N < max_threads) {
        int n_groups = N;
        for (int i = 0; i < n_groups; ++i) {
            groups.push_back({i});
        }
    } else {
        // Distribute the levels uniformly between all cores.
        groups = std::vector<std::vector<int>>(max_threads);
        int i = 0;
        int k = 0;
        while (k < N) {
            groups[i].push_back(k);
            ++k;
            if (i == max_threads - 1) {
                i = 0;
            } else {
                ++i;
            }
        }
    }

    auto perform_warpings = [](auto& source_peaks, auto source_rt_start,
                               auto source_rt_end, auto ref_rt_start,
                               auto ref_rt_end) {
        auto warped_peaks =
            peaks_in_rt_range(source_peaks, source_rt_start, source_rt_end);

        // Warp the peaks by linearly interpolating their retention time
        // to the current segment's. Note that we are just performing
        // linear displacement of the center of the peaks, we do not
        // deform the peak shape by adjusting the sigmas.
        for (auto& peak : warped_peaks) {
            double x =
                (peak.rt - source_rt_start) / (source_rt_end - source_rt_start);
            peak.rt = lerp(ref_rt_start, ref_rt_end, x);
        }
        return warped_peaks;
    };

    auto compute_similarities = [&perform_warpings](
                                    auto& level, auto rt_start, auto rt_end,
                                    auto rt_min, auto delta_rt,
                                    auto& target_peaks, auto& source_peaks) {
        auto target_peaks_segment =
            peaks_in_rt_range(target_peaks, rt_start, rt_end);

        for (auto& warping : level.potential_warpings) {
            int x_start = warping.x_start;
            int x_end = warping.x_end;

            double sample_rt_start = rt_min + warping.x_start * delta_rt;
            double sample_rt_width = (x_end - x_start) * delta_rt;
            double sample_rt_end = sample_rt_start + sample_rt_width;

            auto source_peaks_warped = perform_warpings(
                source_peaks, sample_rt_start, sample_rt_end, rt_start, rt_end);

            // Calculate the peak overlap between the reference and warped
            // peaks for this segment.
            double similarity = Warp2D::similarity_2D(target_peaks_segment,
                                                      source_peaks_warped);
            warping.warped_similarity = similarity;
        }
    };

    std::cout << "finding similarities..." << std::endl;
    std::vector<std::thread> threads(groups.size());
    for (size_t i = 0; i < groups.size(); ++i) {
        threads[i] =
            std::thread([i, &groups, &levels, rt_min, delta_rt,
                         segment_rt_width, &target_peaks_filtered,
                         &source_peaks_filtered, &compute_similarities]() {
                for (const auto& k : groups[i]) {
                    auto& current_level = levels[k];
                    double rt_start = rt_min + k * segment_rt_width;
                    double rt_end = rt_start + segment_rt_width;
                    compute_similarities(
                        current_level, rt_start, rt_end, rt_min, delta_rt,
                        target_peaks_filtered, source_peaks_filtered);
                }
            });
    }
    for (auto& thread : threads) {
        thread.join();
    }

    std::cout << "walking back nodes..." << std::endl;
    for (int k = N - 1; k >= 0; --k) {
        auto& current_level = levels[k];
        const auto& next_level = levels[k + 1];

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

    // Walk back nodes to find optimal warping path.
    std::cout << "finding optimal path..." << std::endl;
    std::vector<int> warp_by;
    warp_by.reserve(N + 1);
    warp_by.push_back(0);
    for (int i = 0; i < N; ++i) {
        auto u = levels[i].nodes[warp_by[i]].u;
        warp_by.push_back(u);
    }

    std::cout << "warping optimal nodes..." << std::endl;

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

        auto warped_peaks_segment = perform_warpings(
            source_peaks, sample_rt_start, sample_rt_end, rt_start, rt_end);
        for (const auto& peak : warped_peaks_segment) {
            warped_peaks.push_back(peak);
        }
    }

    return warped_peaks;
}
