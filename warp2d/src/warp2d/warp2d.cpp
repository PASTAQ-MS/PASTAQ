#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

#include "warp2d/warp2d.hpp"

// TODO(alex): move to header?
struct Node {
    double f;  // Cummulative similarity value.
    int u;     // Optimum predecessor position.
};

// TODO(alex): move to header?
struct Level {
    int start;
    int end;
    int x;
    std::vector<Node> nodes;
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

double Warp2D::peak_overlap(const Centroid::Peak& peak_a,
                            const Centroid::Peak& peak_b) {
    // TODO(alex): Use rt/mz/height or rt_centroid/mz_centroid/height_centroid?

    // Early return if the peaks do not intersect in the +/-2 * sigma_mz/rt
    {
        double min_rt_a = peak_a.rt - 2 * peak_a.sigma_rt;
        double max_rt_a = peak_a.rt + 2 * peak_a.sigma_rt;
        double min_mz_a = peak_a.mz - 2 * peak_a.sigma_mz;
        double max_mz_a = peak_a.mz + 2 * peak_a.sigma_mz;
        double min_rt_b = peak_b.rt - 2 * peak_b.sigma_rt;
        double max_rt_b = peak_b.rt + 2 * peak_b.sigma_rt;
        double min_mz_b = peak_b.mz - 2 * peak_b.sigma_mz;
        double max_mz_b = peak_b.mz + 2 * peak_b.sigma_mz;

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
        filtered_peaks.push_back(peak);
    }
    return filtered_peaks;
}

// TODO(alex): Make sure to verify the integer sizes used here as well as the
// discrepancy between using size_t and int.
std::vector<Centroid::Peak> Warp2D::warp_peaks(
    const std::vector<Centroid::Peak>& target_peaks,
    const std::vector<Centroid::Peak>& source_peaks,
    const Warp2D::Parameters& parameters) {
    // TODO(alex): Use only abbreviated forms in this function (Below).
    int sample_length = parameters.num_points;
    int segment_length = parameters.window_size;
    int slack = parameters.slack;

    // The number of segments.
    int num_segments = sample_length / segment_length;
    sample_length = num_segments * segment_length;

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
    // TODO(alex): Verify that the rt_min/max range and delta_rt correspond with
    // the sample length, etc.
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
    // rt_min = -681;
    // rt_max = 8762.92;
    double delta_rt = (rt_max - rt_min) / (double)(sample_length);
    double rt_sample_width = delta_rt * segment_length;
    // Adjust rt_max to fit all segments.
    rt_max = rt_min + rt_sample_width * num_segments;

    // Filter the peaks in each segment.
    int n_peaks_per_segment = 50;  // FIXME: Hardcoding this for now.
    std::vector<Centroid::Peak> target_peaks_filtered;
    std::vector<Centroid::Peak> source_peaks_filtered;
    for (int i = 0; i < num_segments; ++i) {
        double rt_start = rt_min + i * rt_sample_width;
        double rt_end = rt_start + rt_sample_width;
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
    std::vector<Level> levels(num_segments + 1);
    int N = num_segments;
    int t = slack;
    int m = segment_length;
    int Lt = sample_length;
    for (int i = 0; i < num_segments; ++i) {
        int start = std::max((i * (m - t)), (Lt - (N - i) * (m + t)));
        int end = std::min((i * (m + t)), (Lt - (N - i) * (m - t)));
        int length = end - start + 1;
        levels[i].start = start;
        levels[i].end = end;
        levels[i].nodes = std::vector<Node>(length);
        levels[i].x = i * m;
        for (int j = 0; j < length; ++j) {
            levels[i].nodes[j].f = -std::numeric_limits<double>::infinity();
            levels[i].nodes[j].u = 0;
        }
    }
    levels[num_segments].start = Lt;
    levels[num_segments].end = Lt;
    levels[num_segments].x = Lt;
    levels[num_segments].nodes.push_back({0.0, 0});
    std::cout << "rt_min: " << rt_min << std::endl;
    std::cout << "rt_max: " << rt_max << std::endl;

    // DEBUG
    // int level_i = 0;
    // for (const auto& level : levels) {
    // std::cout << "level: " << level_i << std::endl;
    // std::cout << "start: " << level.start << std::endl;
    // std::cout << "x: " << level.x << std::endl;
    // std::cout << "end: " << level.end << std::endl;
    // std::cout << "n_nodes: " << level.nodes.size() << std::endl;
    // std::cout << std::endl;
    //++level_i;
    //}
    // std::exit(-1);

    // Perform dynamic programming to find optimal warping path.
    for (int i = num_segments - 1; i >= 0; --i) {
        auto& current_level = levels[i];
        const auto& next_level = levels[i + 1];

        double rt_start = rt_min + i * rt_sample_width;
        double rt_end = rt_start + rt_sample_width;
        auto target_peaks_segment =
            peaks_in_rt_range(target_peaks_filtered, rt_start, rt_end);
        auto source_peaks_segment =
            peaks_in_rt_range(source_peaks_filtered, rt_start, rt_end);

        // DEBUG
        // std::cout << "rt_start: " << rt_start << std::endl;
        // std::cout << "rt_end: " << rt_end << std::endl;
        // std::cout << "current_level.start: " << current_level.start
        //<< std::endl;
        // std::cout << "current_level.end: " << current_level.end << std::endl;
        // std::cout << "current_level.x: " << current_level.x << std::endl;
        // std::cout << std::endl;
        std::cout << "level: " << i << " refTimeStart: " << rt_start
                  << " refTimeEnd: " << rt_end << std::endl;

        for (int k = 0; k < (int)current_level.nodes.size(); ++k) {
            int x_i = current_level.start + k;
            auto& node_i = current_level.nodes[k];
            int x_j_min = std::max(x_i + m - t, next_level.start);
            int x_j_max = std::min(x_i + m + t, next_level.end);
            int j_min = x_j_min - next_level.start;
            int j_max = x_j_max - next_level.start;
            for (int j = j_min; j <= j_max; ++j) {
                int x_j = next_level.start + j;
                // std::cout << "x_i: " << x_i << " x_j: " << x_j << std::endl;
                auto& node_j = next_level.nodes[j];

                double warped_time_start = x_i * delta_rt + rt_min;
                double warped_time_end = x_j * delta_rt + rt_min;

                // Make a copy of the peaks for warping.
                std::vector<Centroid::Peak> source_peaks_warped;
                source_peaks_warped.reserve(source_peaks_segment.size());
                for (const auto& peak : source_peaks_segment) {
                    source_peaks_warped.push_back(peak);
                }
                // std::cout << "warped_start: " << warped_start << std::endl;
                // std::cout << "warped_end: " << warped_end << std::endl;

                // Warp the peaks.
                for (auto& peak : source_peaks_warped) {
                    // TODO: Fix to be numerically stable.
                    double lerp_peak =
                        (peak.rt - rt_start) / (rt_end - rt_start) *
                            (warped_time_end - warped_time_start) +
                        warped_time_start;
                    peak.rt = lerp_peak;
                    peak.rt_centroid = lerp_peak;
                }
                // Calculate the peak overlap between the reference and warped
                // peaks.
                double similarity = Warp2D::similarity_2D(target_peaks_segment,
                                                          source_peaks_warped);
                double f_sum = node_j.f + similarity;
                if (f_sum > node_i.f) {
                    node_i.f = f_sum;
                    node_i.u = j;
                }
            }
        }
    }

    // DEBUG
    // int k = 0;
    // for (const auto& level : levels) {
    // std::cout << "level: " << k;
    // for (const auto& node : level.nodes) {
    // std::cout << " [f: " << node.f << " u: " << node.u << "]";
    //}
    //++k;
    // std::cout << std::endl;
    //}

    // Walk back nodes to find optimal warping path.
    std::vector<int> warp_by;
    warp_by.reserve(num_segments + 1);
    warp_by.push_back(0);
    for (int i = 0; i < num_segments; ++i) {
        auto u = levels[i].nodes[warp_by[i]].u;
        warp_by.push_back(u);
        // std::cout << "level: " << i << " u: " << u << std::endl;
    }

    // Warp the sample peaks based on the optimal path.
    std::vector<Centroid::Peak> warped_peaks;
    warped_peaks.reserve(source_peaks.size());
    for (int i = 0; i < num_segments; ++i) {
        double rt_start = rt_min + i * rt_sample_width;
        double rt_end = rt_start + rt_sample_width;
        auto source_peaks_segment =
            peaks_in_rt_range(source_peaks, rt_start, rt_end);

        int x_i = warp_by[i] + levels[i].start;
        int x_j = warp_by[i + 1] + levels[i + 1].start;
        double warped_time_start = x_i * delta_rt + rt_min;
        double warped_time_end = x_j * delta_rt + rt_min;

        // std::cout << "warped_start: " << warped_start << std::endl;
        // std::cout << "warped_end: " << warped_end << std::endl;

        // Warp the peaks.
        for (auto& peak : source_peaks_segment) {
            // TODO: Fix to be numerically stable.
            double lerp_peak = (peak.rt - rt_start) / (rt_end - rt_start) *
                                   (warped_time_end - warped_time_start) +
                               warped_time_start;
            peak.rt = lerp_peak;
            peak.rt_centroid = lerp_peak;
            warped_peaks.push_back(peak);
        }
        // std::cout << "level: " << i << " x_i: " << x_i << " x_j: " << x_j <<
        // std::endl; std::cout << "level: " << i << " levels[i].start: " <<
        // levels[i].start; std::cout << "level: " << i << " levels[i + 1].start:
        // " << levels[i + 1].start;

        // auto u = warp_by[i];
        // auto u_next = warp_by[i + 1];

        ////// Warp the peaks.
        //// double time_diff = u * delta_rt;
        //// for (auto& peak : source_peaks_segment) {
        //// peak.rt += time_diff;
        //// peak.rt_centroid += time_diff;
        //// warped_peaks.push_back(peak);
        ////}

        // double warped_time_end =
        // rt_end + (u_next * delta_rt);  // + optimal warping u?
        //// NOTE: Is this what we want instead:
        //// double warped_time_end = rt_end + (previous_level_u *
        //// delta_rt);
        // double warped_time_start = rt_start + (u * delta_rt);
        //// std::cout << "rt_start: " << rt_start << std::endl;
        //// std::cout << "rt_end: " << rt_end << std::endl;
        //// std::cout << "u: " << u << std::endl;
        //// std::cout << "u_next: " << u_next << std::endl;
        //// std::cout << "warped_start: " << warped_time_start << std::endl;
        //// std::cout << "warped_end: " << warped_time_end << std::endl;
        //// std::cout << "time map: " << warped_time_start << "  " << rt_start
        ////<< std::endl;

        //// Warp the peaks.
        // double time_diff = u * delta_rt;
        // for (auto& peak : source_peaks_segment) {
        //// std::cout << "current rt: " << peak.rt << std::endl;
        //// std::cout << "constant displacement: " << time_diff
        ////<< std::endl;
        //// TODO: Fix to be numerically stable.
        //// double lerp_peak = (peak.rt - rt_start) / (rt_end - rt_start) *
        ////(warped_time_end - warped_time_start) +
        //// warped_time_start;
        // double lerp_peak = (peak.rt - rt_start) / (rt_end - rt_start) *
        //(warped_time_end - warped_time_start) +
        // warped_time_start;
        //// std::cout << "cd_peak: " << peak.rt + time_diff
        ////<< std::endl;
        //// std::cout << "lerp_peak: " << lerp_peak << std::endl;
        //// peak.rt += time_diff;
        //// peak.rt_centroid += time_diff;
        // peak.rt = lerp_peak;
        // peak.rt_centroid = lerp_peak;
        // warped_peaks.push_back(peak);
        //}
    }

    // DEBUG
    std::cout << "DEBUG:" << std::endl;
    auto target_peaks_copy = target_peaks;
    auto source_peaks_copy = source_peaks;
    auto warped_peaks_copy = warped_peaks;
    auto a = filter_peaks(target_peaks_copy, 1000);
    auto b = filter_peaks(source_peaks_copy, 1000);
    auto c = filter_peaks(warped_peaks_copy, 1000);
    double similarity_baseline = Warp2D::similarity_2D(a, a);
    double similarity_before = Warp2D::similarity_2D(a, b);
    double similarity_after = Warp2D::similarity_2D(a, c);
    std::cout << "Similarity Baseline: " << similarity_baseline << std::endl;
    std::cout << "Similarity Before: " << similarity_before << std::endl;
    std::cout << "Similarity After: " << similarity_after << std::endl;
    std::cout << "Similarity Before (Norm): "
              << similarity_before / similarity_baseline << std::endl;
    std::cout << "Similarity After (Norm): "
              << similarity_after / similarity_baseline << std::endl;
    std::cout << "------------" << std::endl;

    return warped_peaks;
}
