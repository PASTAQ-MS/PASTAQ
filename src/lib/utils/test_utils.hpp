#include <vector>
#include <algorithm>  // <-- Needed for std::max_element
#include <numeric>
#include <random>

#include "raw_data/raw_data.hpp"
#include "utils/search.hpp"

using namespace RawData;


RawData::RawData create_test_data(bool block_mode = false) {
    RawData::RawData data;

    // Parameters
    data.instrument_type = Instrument::Type::ORBITRAP;
    data.min_mz = 600.0;
    data.max_mz = 606.0;
    data.min_rt = 0.0;
    data.max_rt = 4.0;
    data.reference_mz = 200;
    data.resolution_ms1 = 17000;
    data.resolution_msn = 17500;
    data.fwhm_rt = 5.0;

    int num_scans = 5;

    // m/z grid parameters
    double mz_start = 602.0;
    double mz_end = 606.0;
    double mz_step = block_mode ? 0.02 : 1.0;  // high-resolution if block_mode
    int num_mz = static_cast<int>(mz_end - mz_start) + 1;

    // Random number generation setup
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(1000.0, 1005.0);

    for (int rt_idx = 0; rt_idx < num_scans; ++rt_idx) {
        Scan scan;
        scan.scan_number = rt_idx;
        scan.ms_level = 1;
        scan.retention_time = static_cast<double>(rt_idx);
        scan.polarity = Polarity::Type::BOTH;

        for (int mz_idx = 0; mz_idx < num_mz; ++mz_idx) {
            double mz_val = mz_start + mz_idx * mz_step;
            scan.mz.push_back(mz_val);

            double intensity = dist(gen);

            if (!block_mode) {
                // Vertical line around m/z = 603
                if (!(rt_idx >= 1 && rt_idx <= 3 && mz_idx == 3)) {
                    intensity = 0.0;
                }
            } else {
                // 3x3 block centered around (603, RT=2)
                if (!(rt_idx >= 1 && rt_idx <=3 && mz_idx >= 2 && mz_idx <= 4)) {
                    intensity = 0.0;
                }
            }

            scan.intensity.push_back(intensity);
        }

        scan.num_points = scan.mz.size();
        scan.max_intensity = *std::max_element(scan.intensity.begin(), scan.intensity.end());
        scan.total_intensity = std::accumulate(scan.intensity.begin(), scan.intensity.end(), 0.0);

        data.scans.push_back(scan);
        data.retention_times.push_back(scan.retention_time);
    }

    return data;
}

RawData::RawData create_test_data_pairs(double initial_separation) {
    RawData::RawData data;

    // Parameters
    data.instrument_type = Instrument::Type::ORBITRAP;
    data.min_mz = 600.0;
    data.max_mz = 606.0;
    data.min_rt = 0.0;
    data.max_rt = 4.0;
    data.reference_mz = 200;
    data.resolution_ms1 = 17000;
    data.resolution_msn = 17500;
    data.fwhm_rt = 5.0;

    int num_scans = 5;
    double mz_step = 0.005;  // fine grid step!

    int num_mz = static_cast<int>((data.max_mz - data.min_mz) / mz_step) + 1;

    // Random number generation
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(1000.0, 1005.0);

    for (int rt_idx = 0; rt_idx < num_scans; ++rt_idx) {
        Scan scan;
        scan.scan_number = rt_idx;
        scan.ms_level = 1;
        scan.retention_time = static_cast<double>(rt_idx);
        scan.polarity = Polarity::Type::BOTH;

        // Generate fine m/z grid
        for (int mz_idx = 0; mz_idx < num_mz; ++mz_idx) {
            double mz_val = data.min_mz + mz_idx * mz_step;
            scan.mz.push_back(mz_val);
            scan.intensity.push_back(0.0);  // initialize with zeros
        }

        double center = 601.0;
        double separation = initial_separation;

        while (center < 606.0 && separation >= 0.01) {
            double pos1 = center - separation / 2.0;
            double pos2 = center + separation / 2.0;

            // Find closest m/z bins and set intensities
            for (size_t idx = 0; idx < scan.mz.size(); ++idx) {
                double mz = scan.mz[idx];
                if (std::abs(mz - pos1) < mz_step/2.0 || std::abs(mz - pos2) < mz_step/2.0) {
                    scan.intensity[idx] = dist(gen);
                }
            }

            center += 1.0;
            separation /= 2.0;
        }

        scan.num_points = scan.mz.size();
        scan.max_intensity = *std::max_element(scan.intensity.begin(), scan.intensity.end());
        scan.total_intensity = std::accumulate(scan.intensity.begin(), scan.intensity.end(), 0.0);

        data.scans.push_back(scan);
        data.retention_times.push_back(scan.retention_time);
    }

    return data;
}

RawData::RawData create_test_data_triplets(double initial_separation) {
    RawData::RawData data;

    // Parameters
    data.instrument_type = Instrument::Type::ORBITRAP;
    data.min_mz = 600.0;
    data.max_mz = 606.0;
    data.min_rt = 0.0;
    data.max_rt = 4.0;
    data.reference_mz = 200;
    data.resolution_ms1 = 17000;
    data.resolution_msn = 17500;
    data.fwhm_rt = 5.0;

    int num_scans = 5;
    double mz_step = 0.005;  // fine m/z grid

    int num_mz = static_cast<int>((data.max_mz - data.min_mz) / mz_step) + 1;

    // Random number generation
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dist(1000.0, 1005.0);

    for (int rt_idx = 0; rt_idx < num_scans; ++rt_idx) {
        Scan scan;
        scan.scan_number = rt_idx;
        scan.ms_level = 1;
        scan.retention_time = static_cast<double>(rt_idx);
        scan.polarity = Polarity::Type::BOTH;

        // Generate fine m/z grid
        for (int mz_idx = 0; mz_idx < num_mz; ++mz_idx) {
            double mz_val = data.min_mz + mz_idx * mz_step;
            scan.mz.push_back(mz_val);
            scan.intensity.push_back(0.0);  // initialize intensities with zeros
        }

        double center = 601.0;
        double separation = initial_separation;

        while (center < 606.0 && separation >= 0.01) {
            double delta = separation / 2.0;

            double pos1 = center - delta;
            double pos2 = center;           // exact center
            double pos3 = center + delta;

            // Find closest m/z bins and set intensities
            for (size_t idx = 0; idx < scan.mz.size(); ++idx) {
                double mz = scan.mz[idx];
                if (std::abs(mz - pos1) < mz_step/2.0 || 
                    std::abs(mz - pos2) < mz_step/2.0 ||
                    std::abs(mz - pos3) < mz_step/2.0) {
                    scan.intensity[idx] = dist(gen);
                }
            }

            center += 1.0;
            separation /= 2.0;
        }

        scan.num_points = scan.mz.size();
        scan.max_intensity = *std::max_element(scan.intensity.begin(), scan.intensity.end());
        scan.total_intensity = std::accumulate(scan.intensity.begin(), scan.intensity.end(), 0.0);

        data.scans.push_back(scan);
        data.retention_times.push_back(scan.retention_time);
    }

    return data;
}

// RawData::RawData create_test_data(bool block_mode = false) {
//     RawData::RawData data;

//     // Parameters
//     data.instrument_type = Instrument::Type::ORBITRAP;
//     data.min_mz = 600.0;
//     data.max_mz = 606.0;
//     data.min_rt = 0.0;
//     data.max_rt = 4.0;
//     data.reference_mz = 200;
//     data.resolution_ms1 = 17000;
//     data.resolution_msn = 17500;
//     data.fwhm_rt = 5.0;

//     int num_scans = 5;
//     int num_mz = 7;

//     // Random number generation setup
//     std::random_device rd;
//     std::mt19937 gen(rd());
//     std::uniform_real_distribution<> dist(1000.0, 1005.0);

//     for (int rt_idx = 0; rt_idx < num_scans; ++rt_idx) {
//         Scan scan;
//         scan.scan_number = rt_idx;
//         scan.ms_level = 1;
//         scan.retention_time = static_cast<double>(rt_idx);
//         scan.polarity = Polarity::Type::BOTH;

//         for (int mz_idx = 600; mz_idx <= (600 + num_mz); ++mz_idx) {
//             double mz_val = static_cast<double>(mz_idx);
//             // double mz_val = 602.5;
//             scan.mz.push_back(static_cast<double>(mz_val));

//             double intensity = dist(gen);

//             // If block_mode is false and this point is not in the vertical line, make it zero
//             if (!block_mode && !(rt_idx >= 1 && rt_idx <= 3 && mz_idx == 603)) {
//                 intensity = 0.0;
//             }

//             // If block_mode is true and this point is not in the 3x3 block, make it zero
//             if (block_mode && !(rt_idx >= 1 && rt_idx <= 3 && mz_idx >= 602 && mz_idx <= 604)) {
//                 intensity = 0.0;
//             }

//             scan.intensity.push_back(intensity);
//         }

//         scan.num_points = scan.mz.size();
//         scan.max_intensity = *std::max_element(scan.intensity.begin(), scan.intensity.end());
//         scan.total_intensity = std::accumulate(scan.intensity.begin(), scan.intensity.end(), 0.0);

//         data.scans.push_back(scan);
//         data.retention_times.push_back(scan.retention_time);
//     }

//     return data;
// }


// RawData::RawData create_test_data(bool block_mode = false) {
//     RawData::RawData data;

//     // Parameters
//     data.min_mz = 0.0;
//     data.max_mz = 6.0;
//     data.min_rt = 0.0;
//     data.max_rt = 4.0;
//     data.reference_mz = 200;
//     data.resolution_ms1 = 17000;
//     data.resolution_msn = 17500;
//     data.fwhm_rt = 5.0;

//     int num_scans = 5;
//     int num_mz = 7;

//     for (int rt_idx = 0; rt_idx < num_scans; ++rt_idx) {
//         Scan scan;
//         scan.scan_number = rt_idx;
//         scan.ms_level = 1;
//         scan.retention_time = static_cast<double>(rt_idx);
//         scan.polarity = Polarity::Type::POSITIVE;

//         for (int mz_idx = 0; mz_idx < num_mz; ++mz_idx) {
//             scan.mz.push_back(static_cast<double>(mz_idx));

//             double intensity = 0.0;

//             if (block_mode) {
//                 // Central 3x3 block
//                 if (rt_idx >= 1 && rt_idx <= 3 && mz_idx >= 2 && mz_idx <= 4) {
//                     intensity = 2.0;
//                 }
//             } else {
//                 // Vertical line at mz index 3
//                 if (rt_idx >= 1 && rt_idx <= 3 && mz_idx == 3) {
//                     intensity = 2.0;
//                 }
//             }

//             scan.intensity.push_back(intensity);
//         }

//         scan.num_points = scan.mz.size();
//         scan.max_intensity = *std::max_element(scan.intensity.begin(), scan.intensity.end());
//         scan.total_intensity = std::accumulate(scan.intensity.begin(), scan.intensity.end(), 0.0);

//         data.scans.push_back(scan);
//         data.retention_times.push_back(scan.retention_time);
//     }

//     return data;
// }
