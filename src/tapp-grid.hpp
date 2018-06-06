#ifndef TAPP_GRID_GRID_HPP
#define TAPP_GRID_GRID_HPP

#include <optional>

namespace Grid {
// Represent the grid dimensions in index coordinates.
struct Dimensions {
    int n;
    int m;
};

// Represent the real world bounds for this grid.
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
    virtual std::optional<double> at(int i, int j) = 0;

    // Get the real world mass/charge stored in the given index.
    virtual std::optional<double> mz(int i) = 0;

    // Get the real world retention time stored in the given index.
    virtual std::optional<double> rt(int j) = 0;

    // Return the dimensions of the Grid in index coordinates:
    //   i <- [0,N], j <- [0,M]
    virtual Dimensions dim() = 0;

    // Return the bounds of the grid.
    virtual Bounds bounds() = 0;
};
}  // namespace Grid
#endif /* TAPP_GRID_GRID_HPP */
