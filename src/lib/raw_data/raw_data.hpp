#ifndef RAWDATA_RAWDATA_HPP
#define RAWDATA_RAWDATA_HPP

#include <cassert>
#include <cmath>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

// The instrument in which the data was acquired.
namespace Instrument {
enum Type : uint8_t { UNKNOWN = 0, QUAD = 1, TOF = 2, FTICR = 3, ORBITRAP = 4 };
}  // namespace Instrument

// The ionization polarity of the analysis. It is possible to configure some
// instruments to acquire data in alternating polarities (+,-,+...). This
// parameter is meant to be used as a filter of the scans that are going to be
// read from the raw data.
namespace Polarity {
enum Type : uint8_t { UNKNOWN = 0, POSITIVE = 1, NEGATIVE = 2, BOTH = 3 };
}  // namespace Polarity

// This describes the fragmentation method for MS/MS spectra.
namespace ActivationMethod {
enum Type : uint8_t { UNKNOWN = 0, CID = 1, HCD = 2 };
}  // namespace ActivationMethod

// The Xic stores the information for an extracted ion chromatogram.
namespace Xic {
// A Xic can be generated using the summation of all points for each scan,
// or the maximum intensity. When used on the entire dataset, the former is
// called the Total Ion Chromatogram (TIC) whereas the later is called base peak
// chromatogram.
enum Method : uint8_t { UNKNOWN = 0, SUM = 1, MAX = 2 };
struct Xic {
    // Resulting data vectors.
    std::vector<double> retention_time;
    std::vector<double> intensity;
    // The parameters used to generate this Xic.
    Method method;
    double min_mz;
    double max_mz;
    double min_rt;
    double max_rt;
};
}  // namespace Xic

// In this namespace we have access to the data structures for working with raw
// data.
namespace RawData {
// For an MS/MS scan, some information about the precursor scan is saved by the
// instrument.
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
    ActivationMethod::Type activation_method;
    // The total isolation window selected for fragmentation in m/z units.
    double window_wideness;
};

struct Scan {
    // Index for this scan.
    uint64_t scan_number;
    // Type of ms_level of this scan (i.e. MS1/MS2/MSn).
    uint64_t ms_level;
    // How many mz-intensity pairs are contained in this scan.
    uint64_t num_points;
    // Retention time in seconds of this scan;
    double retention_time;
    // This is the actual data of the scan, split into two vectors, mz and
    // intensity. These vectors should have the same size (num_points).
    std::vector<double> mz;
    std::vector<double> intensity;
    // Store the maximum and total intensity for this scan.
    double max_intensity;
    double total_intensity;
    // The ionization polarity for this scan.
    Polarity::Type polarity;
    // In case this is a MSn scan, the precursor information will be stored
    // here.
    // TODO(alex): We might want to make this a smart pointer, as this field is
    // only useful for MSn scans.
    PrecursorInformation precursor_information;
};

// Main structure that hold information about a RawData file.
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
    std::vector<Scan> scans;
    // This information is saved for quick search.
    // TODO: Note that this is unnecessary if our search function is able to
    // search through the `scans` array.
    std::vector<double> retention_times;
};

// Raw data points in a struct of arrays format.
struct RawPoints {
    uint64_t num_points;
    std::vector<double> rt;
    std::vector<double> mz;
    std::vector<double> intensity;
};

// Calculate the extracted ion chromatogram for ROI described by the
// min/max_mz/rt on the given raw_data.
Xic::Xic xic(const RawData &raw_data, double min_mz, double max_mz,
             double min_rt, double max_rt, Xic::Method method);

// Calculate the theoretical FWHM of the peak for the given mz.
double theoretical_fwhm(const RawData &raw_data, double mz);

// Transform the FWHM to sigma assuming a Gaussian distribution.
double fwhm_to_sigma(double fwhm);

// Find the raw data points within the square region defined by min/max_mz/rt.
RawPoints raw_points(const RawData &raw_data, double min_mz, double max_mz,
                     double min_rt, double max_rt);
}  // namespace RawData

// In this namespace we have access to the data structures for working with
// identification data in mzIdentML format.
namespace IdentData {
// FIXME: A lot more documentation is necessary here.
struct SpectrumId {
    std::string id;
    bool pass_threshold;
    bool modifications;
    std::string sequence;
    std::string peptide_id;
    uint8_t charge_state;
    double theoretical_mz;
    double experimental_mz;
    double retention_time;
    int64_t rank;
};

struct CVParam {
    std::string name;
    std::string accession;
    std::string cv_ref;
    std::string value;
};

struct PeptideModification {
    double monoisotopic_mass_delta;
    double average_mass_delta;
    std::string residues;
    int64_t location;
    std::vector<CVParam> cv_params;
};

struct Peptide {
    std::string id;
    std::string sequence;
    std::vector<PeptideModification> modifications;
};

struct DBSequence {
    std::string id;
    std::string value;
};

struct ProteinHypothesis {
    std::string db_sequence_id;
    bool pass_threshold;
    std::vector<std::string> spectrum_ids;
};

struct IdentData {
    std::vector<DBSequence> db_sequences;
    std::vector<Peptide> peptides;
    std::vector<SpectrumId> spectrum_ids;
    std::vector<ProteinHypothesis> protein_hypotheses;
};
}  // namespace IdentData

#endif /* RAWDATA_RAWDATA_HPP */
