#ifndef GRID_GRIDRUNNERS_HPP
#define GRID_GRIDRUNNERS_HPP

#include <vector>

#include "grid.hpp"

namespace Grid::Runners::Serial {

// Perform grid splatting for all given points in serial.
std::vector<double> run(const Grid::Parameters &parameters,
                        const std::vector<Point> &all_points);

}  // namespace Grid::Runners::Serial

namespace Grid::Runners::Parallel {

// Perform grid splatting for all given points in parallel.
std::vector<double> run(uint64_t max_threads, const Parameters &parameters,
                        const std::vector<Point> &all_points);

// Splits the parameters in retention time into a maximum of n_split segments.
std::vector<Parameters> split_segments(const Parameters &original_params,
                                       uint64_t n_splits);

// Assign each of the given points into the appropriate segment.
std::vector<std::vector<Point>> assign_points(
    const std::vector<Parameters> &all_parameters,
    const std::vector<Point> &points);

// Merge overlapping segments into a single data vector. Note that the current
// implementation is very memory inefficient, we should avoid copying the
// array of points into the different groups. We can either store the references
// or move the values.
std::vector<double> merge_segments(
    std::vector<Parameters> &parameters_array,
    std::vector<std::vector<double>> &data_array);

}  // namespace Grid::Runners::Parallel

#endif /* GRID_GRIDRUNNERS_HPP */
