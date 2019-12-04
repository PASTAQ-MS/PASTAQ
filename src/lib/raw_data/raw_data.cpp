#include "raw_data/raw_data.hpp"
#include "utils/search.hpp"

double RawData::theoretical_fwhm(const RawData &raw_data, double mz) {
    double e = 0;
    switch (raw_data.instrument_type) {
        case Instrument::ORBITRAP: {
            e = 1.5;
        } break;
        case Instrument::FTICR: {
            e = 2;
        } break;
        case Instrument::TOF: {
            e = 1;
        } break;
        case Instrument::QUAD: {
            e = 0;
        } break;
        case Instrument::UNKNOWN: {
            assert(false);  // Can't handle unknown instruments.
        } break;
    }
    double mz_ref = raw_data.reference_mz;
    double fwhm_ref = mz_ref / raw_data.resolution_ms1;
    return fwhm_ref * std::pow(mz / mz_ref, e);
}

double RawData::fwhm_to_sigma(double fwhm) {
    return fwhm / (2 * std::sqrt(2 * std::log(2)));
}

Xic::Xic RawData::xic(const RawData &raw_data, double min_mz, double max_mz,
                      double min_rt, double max_rt, Xic::Method method) {
    Xic::Xic result = {};
    result.min_mz = min_mz;
    result.max_mz = max_mz;
    result.min_rt = min_rt;
    result.max_rt = max_rt;
    result.method = method;
    const auto &scans = raw_data.scans;
    if (scans.size() == 0) {
        return result;
    }

    // Find scan indices.
    size_t min_j = Search::lower_bound(raw_data.retention_times, min_rt);
    size_t max_j = scans.size();
    if (scans[min_j].retention_time < min_rt) {
        ++min_j;
    }
    for (size_t j = min_j; j < max_j; ++j) {
        const auto &scan = scans[j];
        if (scan.num_points == 0) {
            continue;
        }
        if (scan.retention_time > max_rt) {
            break;
        }

        double internal_min_mz = min_mz;
        double internal_max_mz = max_mz;
        if (internal_min_mz < scan.mz[0]) {
            internal_min_mz = scan.mz[0];
        }
        if (internal_max_mz > scan.mz[scan.num_points - 1]) {
            internal_max_mz = scan.mz[scan.num_points - 1];
        }
        size_t min_i = Search::lower_bound(scan.mz, min_mz);
        size_t max_i = scan.num_points;
        if (scan.mz[min_i] < min_mz) {
            ++min_i;
        }

        double aggregated_intensity = 0;
        switch (method) {
            case Xic::SUM: {
                // Sum all points in the scan.
                for (size_t i = min_i; i < max_i; ++i) {
                    if (scan.mz[i] > max_mz) {
                        break;
                    }
                    aggregated_intensity += scan.intensity[i];
                }
            } break;
            case Xic::MAX: {
                // Find max point in the scan.
                for (size_t i = min_i; i < max_i; ++i) {
                    if (scan.mz[i] > max_mz) {
                        break;
                    }
                    if (scan.intensity[i] > aggregated_intensity) {
                        aggregated_intensity = scan.intensity[i];
                    }
                }
            } break;
            default: {
                result.method = Xic::UNKNOWN;
                return result;
            } break;
        }

        result.retention_time.push_back(scan.retention_time);
        result.intensity.push_back(aggregated_intensity);
    }
    return result;
}

RawData::RawPoints RawData::raw_points(const RawData &raw_data, double min_mz,
                                       double max_mz, double min_rt,
                                       double max_rt) {
    RawPoints raw_points = {};
    const auto &scans = raw_data.scans;
    if (scans.size() == 0) {
        return raw_points;
    }

    size_t min_j = Search::lower_bound(raw_data.retention_times, min_rt);
    size_t max_j = scans.size();
    if (scans[min_j].retention_time < min_rt) {
        ++min_j;
    }

    for (size_t j = min_j; j < max_j; ++j) {
        const auto &scan = scans[j];
        if (scan.retention_time > max_rt) {
            break;
        }
        if (scan.num_points == 0) {
            continue;
        }

        size_t min_i = Search::lower_bound(scan.mz, min_mz);
        size_t max_i = scan.num_points;
        if (scan.mz[min_i] < min_mz) {
            ++min_i;
        }
        bool scan_not_empty = false;
        for (size_t i = min_i; i < max_i; ++i) {
            if (scan.mz[i] > max_mz) {
                break;
            }
            scan_not_empty = true;
            raw_points.rt.push_back(scan.retention_time);
            raw_points.mz.push_back(scan.mz[i]);
            raw_points.intensity.push_back(scan.intensity[i]);
            ++raw_points.num_points;
        }
        if (scan_not_empty) {
            ++raw_points.num_scans;
        }
    }

    return raw_points;
}
