#ifndef GRID_RAWDATA_HPP
#define GRID_RAWDATA_HPP

#include <cassert>
#include <cmath>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

namespace Instrument {

enum Type : uint8_t { QUAD = 0, TOF = 1, FTICR = 2, ORBITRAP = 3, UNKNOWN = 4 };

}  // namespace Instrument

// In this namespace we have access to the data structures for working with the
// read raw data.
namespace RawData {
enum Polarity : uint8_t { POSITIVE = 0, NEGATIVE = 1, BOTH = 2 };
enum ActivationMethod : uint8_t { UNKNOWN = 0, CID = 1, HCD = 2 };

struct PrecursorInformation {
    // Index for the scan that caused the MSn event.
    uint64_t scan_number;
    // Detected charge for the precursor.
    uint8_t charge;
    // Mass to charge of the MSn event.
    double mz;
    // Intensity of the precursor event.
    double intensity;
    // The activation method for the fragmentation of the MSn event.
    ActivationMethod activation_method;
    // The total isolation window selected for fragmentation in m/z units.
    double window_wideness;
};

struct Scan {
    // Index for this scan.
    uint64_t scan_number;
    // Type of ms_level of this scan (i.e. MS1/MS2/MSn).
    uint64_t ms_level;
    // How many mz-intensity pairs are containd in this scan.
    uint64_t num_points;
    // Retention time in seconds of this scan;
    double retention_time;
    // mz-intensity vectors should have the same size (num_points).
    std::vector<double> mz;
    std::vector<double> intensity;
    // The polarity of the ionization for this scan.
    Polarity polarity;
    // In case this is a MSn scan, the precursor information will be stored
    // here.
    PrecursorInformation precursor_information;
};

struct RawData {
    // The instrument type.
    Instrument::Type instrument_type;
    // Min/max mass to charge range (m/z).
    double min_mz;
    double max_mz;
    // Min/max retention time range (seconds).
    double min_rt;
    double max_rt;
    // Resolution of MS1/MSn at the reference m/z. In this case the resolution
    // is defined as:
    //
    //     R = reference_mz/fwhm_at_reference_mz
    //
    double resolution_ms1;
    double resolution_msn;
    double reference_mz;
    // Average full width half maximum of chromatographic peaks.
    double fwhm_rt;

    // Extracted scans.
    std::vector<::RawData::Scan> scans;
    std::vector<double> retention_times;
    std::vector<double> total_ion_chromatogram;
    std::vector<double> base_peak_chromatogram;

    std::tuple<std::vector<double>, std::vector<double>> xic(
        double min_mz, double max_mz, double min_rt, double max_rt,
        std::string method) const;
};

// Calculate the theoretical FWHM of the peak for the given mz.
double theoretical_fwhm(const RawData &raw_data, double mz);

// Transform the FWHM to sigma assuming a Gaussian distribution.
double fwhm_to_sigma(double fwhm);

}  // namespace RawData

#endif /* GRID_RAWDATA_HPP */
