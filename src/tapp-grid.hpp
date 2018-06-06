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
//   rtRange <- [minRt, maxRt]
//   mzRange <- [minMz, maxMz]
struct Bounds {
    double minRt;
    double maxRt;
    double minMz;
    double maxMz;
};

class Interface {
   public:
    // Fetch the data stored at the given indices.
    // TODO: Make this a [][] overloaded operator, along with a way of storing
    // the data on them.
    virtual std::optional<double> at(unsigned int i, unsigned int j) = 0;

    // Get the real world mass/charge stored in the given index.
    virtual std::optional<double> mz(unsigned int i) = 0;

    // Get the real world retention time stored in the given index.
    virtual std::optional<double> rt(unsigned int j) = 0;

    // Return the dimensions of the Grid in index coordinates:
    //   i <- [0,N], j <- [0,M]
    virtual Dimensions dim() = 0;

    // Return the bounds of the grid.
    virtual Bounds bounds() = 0;
};
}  // namespace Grid

#endif /* TAPP_GRID_GRID_HPP */
