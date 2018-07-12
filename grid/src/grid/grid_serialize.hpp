#ifndef GRID_GRIDSERIALIZE_HPP
#define GRID_GRIDSERIALIZE_HPP

#include <iostream>

#include "grid.hpp"

// This namespace groups the functions used to serialize Grid data structures
// into a binary stream.
namespace Grid::Serialize {

// Grid::Parameters
bool read_parameters(std::istream &stream, Grid::Parameters *parameters);
bool write_parameters(std::ostream &stream, const Grid::Parameters &parameters);

// Grid::Point
bool read_point(std::istream &stream, Grid::Point *point);
bool write_point(std::ostream &stream, const Grid::Point &point);

}  // namespace Grid::Serialize

#endif /* GRID_GRIDSERIALIZE_HPP */
