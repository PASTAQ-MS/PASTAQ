#ifndef GRID_GRID_HPP
#define GRID_GRID_HPP

#include <cstdint>
#include <vector>

#include "raw_data/raw_data.hpp"

namespace Grid {

struct Mesh {
    // Represent the grid dimensions in index coordinates:
    //   n: Number of sampling points in mz.
    //   m: Number of sampling points in rt.
    //   k: Number of sampling points per FWHM in mz.
    //   k: Number of sampling points per FWHM in rt.
    uint64_t n;
    uint64_t m;
    uint64_t k;
    uint64_t t;
    // The Mesh data is stored as an array, and the mz and rt corresponding to
    // each bin is memoized for quick indexing when searching.
    // TODO(alex): Not sure that these should be doubles anymore, it uses
    // a significant amount of memory and is storing smoothed data, which
    // already alters the precision of the initial measurements.
    std::vector<double> data;
    std::vector<double> bins_mz;
    std::vector<double> bins_rt;

    // Parameters extracted from the raw data.
    // TODO: Are these necessary to store here? We have memoized the bins in
    // the respective dimensions, so these might not be used at all.
    Instrument::Type instrument_type;
    double reference_mz;
    double fwhm_mz;
    double fwhm_rt;

    // Represent the bounds of this grid in mass-to-charge and retention time
    // coordinates.
    // TODO: Is storing this necessary? We already have this information from
    // the bins_mz/bins_rt array.
    double min_mz;
    double max_mz;
    double min_rt;
    double max_rt;

    // TODO: Methods vs functions? The eternal question, but I feel we should
    // be consistent across all modules.
    uint64_t x_index(double mz);
    uint64_t y_index(double rt);
    double mz_at(uint64_t i);
    double rt_at(uint64_t j);
};

// Applies 2D kernel smoothing. The smoothing is performed in two passes.  First
// the raw data points are mapped into a 2D matrix by splatting them. Sparse
// areas might result in artifacts when the data is noisy, for this reason, the
// data is smoothed again.
//
// Since multiple passes of a Gaussian smoothing is equivalent to a single pass
// with `sigma = sqrt(2) * sigma_pass`, we adjust the sigmas for each pass
// accordingly.
Mesh resample(const RawData::RawData &raw_data, uint64_t num_samples_mz,
              uint64_t num_samples_rt, double smoothing_coef_mz,
              double smoothing_coef_rt);

}  // namespace Grid

#endif /* GRID_GRID_HPP */
