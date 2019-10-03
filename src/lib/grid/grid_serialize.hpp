#ifndef GRID_GRIDSERIALIZE_HPP
#define GRID_GRIDSERIALIZE_HPP

#include <iostream>

#include "grid.hpp"

// This namespace groups the functions used to serialize Grid data structures
// into a binary stream.
namespace Grid::Serialize {

// Grid::Mesh
bool read_mesh(std::istream &stream, Grid::Mesh *mesh);
bool write_mesh(std::ostream &stream, const Grid::Mesh &mesh);

}  // namespace Grid::Serialize

#endif /* GRID_GRIDSERIALIZE_HPP */
