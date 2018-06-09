#ifndef TAPP_GRID_GRID_HPP
#define TAPP_GRID_GRID_HPP

#include <optional>

namespace Grid {

// Represent the grid dimensions in index coordinates:
//   (n == number of columns == Number of points in mz)
//   (m == number of rows == Number of points in rt)
struct Dimensions {
    unsigned int n;
    unsigned int m;
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

class Interface {
   public:
    // Get the value stored at the given position.
    virtual std::optional<double> value_at(unsigned int i, unsigned int j) = 0;

    // Set the value at the given position. Returns the success or failure of
    // the operation.
    virtual bool set_value(unsigned int i, unsigned int j, double value) = 0;

    // Get the real world mass/charge stored in the given index.
    virtual std::optional<double> mz_at(unsigned int i) = 0;

    // Get the real world retention time stored in the given index.
    virtual std::optional<double> rt_at(unsigned int j) = 0;

    // Get the x index of the closest point (rounded down) for a given mz.
    virtual std::optional<unsigned int> x_index(double mz) = 0;

    // Get the y index of the closest point (rounded down) for a given rt.
    virtual std::optional<unsigned int> y_index(double rt) = 0;

    // Scale the sigma for a given mass based on the configuration of the
    // instrument that generated the data.
    virtual double sigma_at_mz(double mz) = 0;

    // Return the dimensions of the Grid in index coordinates:
    //   i <- [0,N], j <- [0,M]
    virtual Dimensions dim() = 0;

    // Return the bounds of the grid.
    virtual Bounds bounds() = 0;
};

// Perform gaussian splatting of the given point into the grid, returns the
// success or failure of the operation.
bool splash(Grid::Interface &grid, double mz, double rt, double value,
            double sigma_rt);

}  // namespace Grid
namespace Instrument {
enum Type { QUAD, TOF, FTICR, ORBITRAP };

struct Config {
    Type type;
    double mz_at_sigma;
    double sigma;
};
}  // namespace Instrument

#endif /* TAPP_GRID_GRID_HPP */
