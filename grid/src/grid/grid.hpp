#ifndef GRID_GRID_HPP
#define GRID_GRID_HPP

#include <vector>
#include <cstdint>

namespace Instrument {

enum Type : unsigned char { QUAD, TOF, FTICR, ORBITRAP };

}  // namespace Instrument

namespace Grid {

// Flags
enum Flags : unsigned char { WARPED_MESH = 0b00000001 };

// Represent the grid dimensions in index coordinates:
//   (n == number of columns == Number of points in mz)
//   (m == number of rows == Number of points in rt)
struct Dimensions {
    uint64_t n;
    uint64_t m;
};

// Represent the real world bounds for this grid. Note that the bounds are
// included, that is the minimum represent the first point in the grid for that
// dimension and the maximum the last.
//   rt_range <- [min_rt, max_rt]
//   mz_range <- [min_mz, max_mz]
struct Bounds {
    double min_rt;
    double max_rt;
    double min_mz;
    double max_mz;
};

// The parameters necessary to perform gaussian smoothing. Note that since the
// number of sampling points is dependent on the sigma at a particular mz, we
// need to pass this data in aggregate. The number of sampling points for the
// retention time should be more consistent across the run, so there is no need
// to correct for it.
struct SmoothingParams {
    double mz;
    double sigma_mz;
    double sigma_rt;
};

// This groups all parameters that can be used for operating with Grid objects.
struct Parameters {
    Dimensions dimensions;
    Bounds bounds;
    SmoothingParams smoothing_params;
    unsigned char instrument_type;
    unsigned char flags;
};

// A Grid::Point represents a single measured value at a given mz and rt.
struct Point {
    double mz;
    double rt;
    double value;
};

// Perform gaussian splatting of the given point into the grid, returns the
// success or failure of the operation.
bool splat(const Point &point, const Parameters &parameters,
           std::vector<double> &data);

// Get the real world mass/charge stored in the given index for the given
// parameters. It doesn't perform any boundary checks.
double mz_at(uint64_t i, const Parameters &parameters);

// Get the real world retention time stored in the given index. It doesn't
// perform any boundary checks.
double rt_at(uint64_t j, const Parameters &parameters);

// Get the x index of the closest point (rounded down) for a given mz. It
// doesn't perform any boundary checks.
uint64_t x_index(double mz, const Parameters &parameters);

// Get the y index of the closest point (rounded down) for a given rt. It
// doesn't perform any boundary checks.
uint64_t y_index(double rt, const Parameters &parameters);

// Get the sigma_mz used for smoothing. In order to maintain the same number
// of sampling points for smoothing across all the mz range of the
// instrument, we need to scale the sigma accordingly.
double sigma_mz(double mz, const Parameters &parameters);

// Get the sigma_rt used for smoothing.
double sigma_rt(const Parameters &parameters);

// Set up the Grid::Dimensions inside the given Grid::Parameters. The dimensions
// of the grid depend on it being a warped mesh or not, the bounds, the sampling
// delta (As defined on the smoothing parameters) and the instrument.
bool calculate_dimensions(Parameters &parameters);

}  // namespace Grid

#endif /* GRID_GRID_HPP */
