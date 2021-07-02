#ifndef GRID_GRID_HPP
#define GRID_GRID_HPP

#include <cstdint>
#include <vector>

#include "raw_data/raw_data.hpp"

namespace Grid {

struct Grid {
    // Represent the grid dimensions in index coordinates:
    //   n: Number of sampling points in mz.
    //   m: Number of sampling points in rt.
    //   k: Number of sampling points per FWHM in mz.
    //   t: Number of sampling points per FWHM in rt.
    uint64_t n;
    uint64_t m;
    uint64_t k;
    uint64_t t;

    // The Grid data is stored as an array, and the mz and rt corresponding to
    // each bin is memoized for quick indexing when searching.
    // TODO(alex): Not sure that these should be doubles anymore, it uses
    // a significant amount of memory and is storing smoothed data, which
    // already alters the precision of the initial measurements. The issue is
    // that the resampling procedure involves the summation of many numbers on
    // several grid points, which might significantly alter results. More
    // testing is needed.
    std::vector<double> data;
    std::vector<double> bins_mz;
    std::vector<double> bins_rt;

    // Parameters extracted from the raw data.
    Instrument::Type instrument_type;
    double reference_mz;
    double fwhm_mz;
    double fwhm_rt;

    // Represent the bounds of this grid in mass-to-charge and retention time
    // coordinates.
    double min_mz;
    double max_mz;
    double min_rt;
    double max_rt;
};

// Applies 2D kernel smoothing. The smoothing is performed in two passes.  First
// the raw data points are mapped into a 2D matrix by splatting them. Sparse
// areas might result in artifacts when the data is noisy, for this reason, the
// data is smoothed again.
//
// Since multiple passes of a Gaussian smoothing is equivalent to a single pass
// with `sigma = sqrt(2) * sigma_pass`, we adjust the sigmas for each pass
// accordingly.
struct ResampleParams {
    uint64_t num_samples_mz;
    uint64_t num_samples_rt;
    double smoothing_coef_mz;
    double smoothing_coef_rt;
};
Grid resample(const RawData::RawData &raw_data, const ResampleParams &params);

// Calculate the index i/j for the given mz/rt on the grid. This calculation is
// performed in linear time.
uint64_t x_index(const Grid &grid, double mz);
uint64_t y_index(const Grid &grid, double rt);

// Calculate the mz/rt at index i/j for a given grid. The calculation is
// performed in linear time.
double mz_at(const Grid &grid, uint64_t i);
double rt_at(const Grid &grid, uint64_t j);

// Extract a subset from the grid based on the given constrained dimensions.
Grid subset(Grid grid, double min_mz, double max_mz, double min_rt, double max_rt);

}  // namespace Grid

#endif /* GRID_GRID_HPP */
