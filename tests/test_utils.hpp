#ifndef TESTS_TESTUTILS
#define TESTS_TESTUTILS
#include <math.h>
#include <cstdint>

namespace TestUtils {
// Check the approximate equality between two floating point numbers. The
// results are truncated to the closest integer for the given precision.
inline bool compare_double(double a, double b, uint64_t precision = 4) {
    double exp = std::pow(10.0, precision);
    return (int64_t)(a * exp) == (int64_t)(b * exp);
}
}  // namespace TestUtils

#endif /* TESTS_TESTUTILS */
