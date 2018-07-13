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
        peaks.push_back(peak);
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
        // Perform splatting in this group.
        threads[i] =
            std::thread([&groups, &parameters, &peaks_array, &data, i]() {
                // Building peaks.
                for (const auto &point : groups[i]) {
                    auto peak = Centroid::build_peak(point, parameters, data);
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
    return peaks;
}

// std::vector<double> Grid::Runners::Parallel::merge_segments(
// std::vector<Grid::Parameters> &parameters_array,
// std::vector<std::vector<double>> &data_array) {
// std::vector<double> merged;
//// Early return if there are errors.
// if (data_array.empty() || parameters_array.empty() ||
// parameters_array.size() != data_array.size()) {
// return merged;
//}
// merged.insert(end(merged), begin(data_array[0]), end(data_array[0]));

// for (size_t n = 1; n < data_array.size(); ++n) {
// auto beg_next = 1 + Grid::y_index(parameters_array[n - 1].bounds.max_rt,
// parameters_array[n]);
//// Sum the overlapping sections.
// int i = merged.size() - beg_next * parameters_array[n - 1].dimensions.n;
// for (size_t j = 0;
// j < (parameters_array[n - 1].dimensions.n * beg_next); ++j) {
// merged[i] += data_array[n][j];
//++i;
//}
//// Insert the next slice.
// merged.insert(
// end(merged),
// begin(data_array[n]) + beg_next * parameters_array[n].dimensions.n,
// end(data_array[n]));
//}
// return merged;
//}

//// std::vector<double> Grid::Runners::Parallel::run(
//// uint64_t max_threads, const Grid::Parameters &parameters,
//// const std::vector<Grid::Point> &all_points) {
////// Split parameters and points into the corresponding  groups.
//// auto all_parameters = split_segments(parameters, max_threads);
//// auto groups = assign_points(all_parameters, all_points);

////// Splatting!
//// std::vector<std::thread> threads(groups.size());
//// std::vector<std::vector<double>> data_array(groups.size());
//// for (size_t i = 0; i < groups.size(); ++i) {
////// Allocate data memory.
//// data_array[i] = std::vector<double>(all_parameters[i].dimensions.n *
//// all_parameters[i].dimensions.m);

////// Perform splatting in this group.
//// threads[i] = std::thread([&groups, &all_parameters, &data_array, i]() {
//// for (const auto &point : groups[i]) {
//// Grid::splat(point, all_parameters[i], data_array[i]);
////}
////});
////}

////// Wait for the threads to finish.
//// for (auto &thread : threads) {
//// thread.join();
////}

//// return merge_segments(all_parameters, data_array);
////}

//// std::vector<Grid::Parameters> Grid::Runners::Parallel::split_segments(
//// const Grid::Parameters &original_params, uint64_t n_splits) {
////// In order to determine the overlapping of the splits we need to
////// calculate what is the maximum distance that will be used by the
////// kernel smoothing. To avoid aliasing we will overlap at least the
////// maximum kernel width.
//////
////// The kernel in rt is always 2 * sigma_rt in both directions.
//// uint64_t kernel_width = Grid::y_index(
//// original_params.bounds.min_rt + 4 * Grid::sigma_rt(original_params),
//// original_params);

////// We need to make sure that we have the minimum number of points for
////// the splits. Since we have an overlap of a single kernel_width, we
////// need to have at least twice that amount of points in order to support
////// full overlap in both directions.
//// uint64_t min_segment_width = 2 * kernel_width;
//// uint64_t segment_width = original_params.dimensions.m / n_splits;
//// if (segment_width < min_segment_width) {
//// segment_width = min_segment_width;
////}

////// If the orginal parameters don't contain the required minimum number
////// of points for segmentation, we can only use one segment.
//// if ((original_params.dimensions.m - 1) < min_segment_width) {
//// return std::vector<Grid::Parameters>({original_params});
////}

////// How many segments do we have with the given segment_width.
//// uint64_t num_segments = original_params.dimensions.m / segment_width;
//// if (original_params.dimensions.m % segment_width) {
////// If we need more segments that the maximum we specify we need to
////// try one less split and adjust the sizes accordingly.
//// if (num_segments + 1 > n_splits) {
//// segment_width = original_params.dimensions.m / (n_splits - 1);
//// num_segments = original_params.dimensions.m / segment_width;
//// if (original_params.dimensions.m % segment_width) {
////++num_segments;
////}
////}
////}

//// std::vector<Grid::Parameters> all_parameters;
//// for (size_t i = 0; i < num_segments; ++i) {
//// uint64_t min_i = 0;
//// uint64_t max_i = 0;
//// if (i == 0) {
//// min_i = 0;
////} else {
//// min_i = segment_width * i - min_segment_width;
////}
//// max_i = segment_width * (i + 1) - 1;
//// if (max_i > original_params.dimensions.m) {
//// max_i = original_params.dimensions.m - 1;
////}
////// Prepare the next Grid::Parameters object.
//// Grid::Parameters parameters(original_params);
//// parameters.bounds.min_rt = Grid::rt_at(min_i, original_params);
//// parameters.bounds.max_rt = Grid::rt_at(max_i, original_params);
//// parameters.dimensions.m = max_i - min_i + 1;
//// all_parameters.emplace_back(parameters);
////}
//// return all_parameters;
////}

//// std::vector<std::vector<Grid::Point>>
/// Grid::Runners::Parallel::assign_points( / const
/// std::vector<Grid::Parameters> &all_parameters, / const
/// std::vector<Grid::Point> &points) { / std::vector<std::vector<Grid::Point>>
/// groups(all_parameters.size());

//// for (const auto &point : points) {
//// for (size_t i = 0; i < all_parameters.size(); ++i) {
//// auto parameters = all_parameters[i];
//// double sigma_rt = Grid::sigma_rt(parameters);
//// if (point.rt + 2 * sigma_rt < parameters.bounds.max_rt) {
//// groups[i].push_back(point);
//// break;
////}
////// If we ran out of parameters this point is assigned to the
////// last one.
//// if (i == all_parameters.size() - 1) {
//// groups[i].push_back(point);
////}
////}
////}
//// return groups;
////}

//// std::vector<double> Grid::Runners::Parallel::merge_segments(
//// std::vector<Grid::Parameters> &parameters_array,
//// std::vector<std::vector<double>> &data_array) {
//// std::vector<double> merged;
////// Early return if there are errors.
//// if (data_array.empty() || parameters_array.empty() ||
//// parameters_array.size() != data_array.size()) {
//// return merged;
////}
//// merged.insert(end(merged), begin(data_array[0]), end(data_array[0]));

//// for (size_t n = 1; n < data_array.size(); ++n) {
//// auto beg_next = 1 + Grid::y_index(parameters_array[n - 1].bounds.max_rt,
//// parameters_array[n]);
////// Sum the overlapping sections.
//// int i = merged.size() - beg_next * parameters_array[n - 1].dimensions.n;
//// for (size_t j = 0;
//// j < (parameters_array[n - 1].dimensions.n * beg_next); ++j) {
//// merged[i] += data_array[n][j];
////++i;
////}
////// Insert the next slice.
//// merged.insert(
//// end(merged),
//// begin(data_array[n]) + beg_next * parameters_array[n].dimensions.n,
//// end(data_array[n]));
////}
//// return merged;
////}
