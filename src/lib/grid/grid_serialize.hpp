#ifndef GRID_GRIDSERIALIZE_HPP
#define GRID_GRIDSERIALIZE_HPP

#include <iostream>

#include "grid.hpp"

// This namespace groups the functions used to serialize Grid data structures
// into a binary stream.
namespace Grid::Serialize {

// Grid::Grid
bool read_grid(std::istream &stream, Grid *grid);
bool write_grid(std::ostream &stream, const Grid &grid);

}  // namespace Grid::Serialize

#endif /* GRID_GRIDSERIALIZE_HPP */
