#include "utils/interpolation.hpp"

double Interpolation::lerp(double y_0, double y_1, double x) {
    return (1 - x) * y_0 + x * y_1;
}
