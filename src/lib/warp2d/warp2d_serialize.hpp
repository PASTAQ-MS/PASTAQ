#ifndef WARP2D_WARP2DSERIALIZE_HPP
#define WARP2D_WARP2DSERIALIZE_HPP

#include <iostream>

#include "warp2d/warp2d.hpp"

// This namespace groups the functions used to serialize Warp2D data
// structures into a binary stream.
namespace Warp2D::Serialize {

bool read_time_map(std::istream &stream, TimeMap *time_map);
bool write_time_map(std::ostream &stream, const TimeMap &time_map);

}  // namespace Warp2D::Serialize

#endif /* WARP2D_WARP2DSERIALIZE_HPP */
