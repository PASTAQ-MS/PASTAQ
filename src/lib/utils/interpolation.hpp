#ifndef UTILS_INTERPOLATION_HPP
#define UTILS_INTERPOLATION_HPP

namespace Interpolation {

// Perform numerically stable linear interpolation of x between y_0 and y_1.
// Note that x is a number between 0 and 1. This is the equivalent of the
// following formula:
//
//     (y - y_0) / (y_1 - y_0) = x;
//     y = x * (y_1 - y_0) + y_0;
//
double lerp(double y_0, double y_1, double x);

}  // namespace Interpolation

#endif /* UTILS_INTERPOLATION_HPP */
