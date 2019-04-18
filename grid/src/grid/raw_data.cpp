#include "grid/raw_data.hpp"

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
