#ifndef TESTS_TESTUTILS
#define TESTS_TESTUTILS
#include <math.h>
#include <cstdint>

#include "centroid/centroid.hpp"

namespace TestUtils {
// Check the approximate equality between two floating point numbers. The
// results are truncated to the closest integer for the given precision.
inline bool compare_double(double a, double b, uint64_t precision = 4) {
    double exp = std::pow(10.0, precision);
    return (int64_t)(a * exp) == (int64_t)(b * exp);
}

inline Centroid::Peak mock_gaussian_peak(size_t id, double height, double mz,
                                         double rt, double sigma_mz,
                                         double sigma_rt) {
    Centroid::Peak peak = {};
    peak.id = id;
    peak.local_max_mz = mz;
    peak.local_max_rt = rt;
    peak.local_max_height = height;
    peak.rt_delta = 0.0;
    peak.roi_min_mz = mz - 3 * sigma_mz;
    peak.roi_max_mz = mz + 3 * sigma_mz;
    peak.roi_min_rt = rt - 3 * sigma_rt;
    peak.roi_max_rt = rt + 3 * sigma_rt;
    peak.raw_roi_mean_mz = mz;
    peak.raw_roi_mean_rt = rt;
    peak.raw_roi_sigma_mz = sigma_mz;
    peak.raw_roi_sigma_rt = sigma_rt;
    peak.raw_roi_max_height = height;
    peak.raw_roi_total_intensity = height;
    peak.fitted_height = height;
    peak.fitted_mz = mz;
    peak.fitted_rt = rt;
    peak.fitted_sigma_mz = sigma_mz;
    peak.fitted_sigma_rt = sigma_rt;
    return peak;
}

}  // namespace TestUtils

#endif /* TESTS_TESTUTILS */
