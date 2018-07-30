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

std::vector<Centroid::Peak> Warp2D::warp_peaks(
    const std::vector<Centroid::Peak>& target_peaks,
    const std::vector<Centroid::Peak>& source_peaks, size_t sample_length,
    size_t segment_length, size_t slack) {
    // The number of segments.
    size_t num_segments = sample_length / segment_length;

    // Minimum time step.
    // DEBUG: Hardcoding values here.
    // double rt_min = 0.0;
    // double rt_max = 40.0;
    // double delta_rt = (rt_max - rt_min) / (double)(sample_length - 1);

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
        for (size_t j = 0; j < length; ++j) {
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
        std::cout << "start: " << current_level.start
                  << " end: " << current_level.end << std::endl;  // DEBUG

        // TODO: Warp peaks here and store in vector.
        // TODO: Calculate similarities here and store in vector.
        std::vector<double> benefit_vector(next_level.nodes.size());

        for (int k = 0; k < current_level.nodes.size(); ++k) {
            auto& node = current_level.nodes[k];
            for (int u = -t; u <= t; ++u) {
                int offset =
                    (current_level.start + k + m + u) - next_level.start;
                if (offset < 0 || offset > next_level.nodes.size() - 1) {
                    continue;
                }

                std::cout << "offset: " << offset << std::endl;
                std::cout << "F(i+1, offset): " << next_level.nodes[offset].f
                          << std::endl;
                std::cout << "U(i+1, offset): " << next_level.nodes[offset].u
                          << std::endl;

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
