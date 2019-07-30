#include "grid/raw_data.hpp"
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

std::tuple<std::vector<double>, std::vector<double>> RawData::RawData::xic(
    double min_mz, double max_mz, double min_rt, double max_rt,
    std::string method) const {
    std::vector<double> rt;
    std::vector<double> intensity;
    const auto &scans = this->scans;
    if (scans.size() == 0) {
        return {rt, intensity};
    }

    // Find scan indices.
    size_t min_j = Search::lower_bound(this->retention_times, min_rt);
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
        if (method == "sum") {
            // Sum all points in the scan.
            for (size_t i = min_i; i < max_i; ++i) {
                if (scan.mz[i] > max_mz) {
                    break;
                }
                aggregated_intensity += scan.intensity[i];
            }
        }
        if (method == "max") {
            // Find max point in the scan.
            for (size_t i = min_i; i < max_i; ++i) {
                if (scan.mz[i] > max_mz) {
                    break;
                }
                if (scan.intensity[i] > aggregated_intensity) {
                    aggregated_intensity = scan.intensity[i];
                }
            }
        }
        // FIXME: Non exhaustive pattern matching. Should we use a string here
        // or an enum?

        rt.push_back(scan.retention_time);
        intensity.push_back(aggregated_intensity);
    }
    return {rt, intensity};
}
