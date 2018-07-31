#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>

#include "warp2d/warp2d.hpp"

struct Node {
    double f;  // Cummulative similarity value.
    int u;     // Optimum predecessor position.
};

struct Level {
    int start;
    int end;
    std::vector<Node> nodes;
};

std::vector<Centroid::Peak> peaks_in_rt_range(
    const std::vector<Centroid::Peak>& source_peaks, double time_start,
    double time_end) {
    std::vector<Centroid::Peak> ret;
    for (const auto& peak : source_peaks) {
        if (peak.rt >= time_start && peak.rt < time_end) {
            ret.push_back(peak);
        }
        // TODO(alex): Should we use rt or rt_centroid?
    }
    return ret;
}

// TODO(alex): Make sure to verify the integer sizes used here as well as the
// discrepancy between using size_t and int.
std::vector<Centroid::Peak> Warp2D::warp_peaks(
    const std::vector<Centroid::Peak>& target_peaks,
    const std::vector<Centroid::Peak>& source_peaks, size_t sample_length,
    size_t segment_length, size_t slack) {
    // The number of segments.
    int num_segments = sample_length / segment_length;

    // Minimum time step.
    // DEBUG: Hardcoding values here.
    // double rt_min = 0.0;
    // double rt_max = 40.0;
    // double delta_rt = (rt_max - rt_min) / (double)(sample_length - 1);

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
    double delta_rt = (rt_max - rt_min) / (double)(sample_length - 1);
    double rt_sample_width = delta_rt * segment_length;
    // Adjust rt_max to fit all segments.
    rt_max = rt_min + rt_sample_width * num_segments;

    // DEBUG
    std::cout << "rt_min: " << rt_min << std::endl;
    std::cout << "rt_max: " << rt_max << std::endl;
    std::cout << "rt_sample_width: " << rt_sample_width << std::endl;
    std::cout << "delta_rt: " << delta_rt << std::endl;

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
        for (int j = 0; j < length; ++j) {
            levels[i].nodes[j].f = -std::numeric_limits<double>::infinity();
            levels[i].nodes[j].u = -1;
        }
    }
    levels[num_segments].start = Lt;
    levels[num_segments].end = Lt;
    levels[num_segments].nodes.push_back({0.0, 0});

    // Perform dynamic programming to find optimal warping path.
    for (int i = num_segments - 1; i >= 0; --i) {
        auto& current_level = levels[i];
        const auto& next_level = levels[i + 1];
        // std::cout << "start: " << current_level.start
        //<< " end: " << current_level.end << std::endl;  // DEBUG
        // TODO: Fetch the peaks belonging to this sector for the reference.
        double rt_start = rt_min + i * rt_sample_width;
        double rt_end = rt_start + rt_sample_width;
        std::cout << "rt_start: " << rt_start << std::endl;
        std::cout << "rt_end: " << rt_end << std::endl;
        auto target_peaks_segment =
            peaks_in_rt_range(target_peaks, rt_start, rt_end);
        auto source_peaks_segment =
            peaks_in_rt_range(source_peaks, rt_start, rt_end);

        // TODO: Warp peaks here and store in vector.
        // TODO: Calculate similarities here and store in vector.
        std::vector<double> benefit_vector(next_level.nodes.size());

        for (int k = 0; k < (int)current_level.nodes.size(); ++k) {
            auto& node = current_level.nodes[k];
            for (int u = -t; u <= t; ++u) {
                int offset = current_level.start - next_level.start + k + m + u;
                if (offset < 0 || offset > (int)next_level.nodes.size() - 1) {
                    continue;
                }

                // DEBUG
                // std::cout << "offset: " << offset << std::endl;
                // std::cout << "F(i+1, offset): " << next_level.nodes[offset].f
                //<< std::endl;
                // std::cout << "U(i+1, offset): " << next_level.nodes[offset].u
                //<< std::endl;

                // Make a copy of the peaks for warping.
                std::vector<Centroid::Peak> source_peaks_copy;
                for (const auto& peak : source_peaks) {
                    source_peaks_copy.push_back(peak);
                }

                // Warp the peaks.
                double time_diff = u * delta_rt;
                for (auto& peak : source_peaks_copy) {
                    peak.rt += time_diff;
                    peak.rt_centroid += time_diff;
                }

                // Calculate the peak overlap between the reference and warped
                // peaks.

                double f_sum =
                    next_level.nodes[offset].f + benefit_vector[offset];
                if (f_sum > node.f) {
                    node.f = f_sum;
                    node.u = u;
                }
            }
        }
    }

    std::vector<Centroid::Peak> warped_peaks;
    return warped_peaks;
}
