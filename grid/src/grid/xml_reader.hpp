#ifndef GRID_XMLREADER_HPP
#define GRID_XMLREADER_HPP

#include <map>
#include <optional>

#include "grid.hpp"

// In this namespace we have access to the data structures for working with the
// read raw data.
namespace RawData {
enum Polarity : unsigned char { POSITIVE, NEGATIVE, BOTH };
enum ActivationMethod : unsigned char { CID, HCD };

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
    PrecursorInformation *precursor_information;
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
};
}  // namespace RawData

// XmlReader handles reading xml tags and the supported xml files. This is not
// a generic module, it is tailored to read Grid::Point objects to use for
// splatting on the grid.
namespace XmlReader {

// An Xml tag.
struct Tag {
    std::string name;
    std::map<std::string, std::string> attributes;
    bool closed;
};

// Reads the next scan that fit into the given parameters from a mzXML file into
// a vector of Grid::Point.
std::optional<std::vector<Grid::Point>> read_next_scan(
    std::istream &stream, const Grid::Parameters &parameters);

// Reads the contents of the next tag in the stream.
std::optional<Tag> read_tag(std::istream &stream);

// Read data until the next tag is found and trim whitespace at the beginning in
// necessary.
std::optional<std::string> read_data(std::istream &stream);

// Read an entire mzxml file into the RawData::RawData data structure filtering
// based on min/max mz/rt and polarity.
std::optional<RawData::RawData> read_mzxml(
    std::istream &stream, double min_mz, double max_mz, double min_rt,
    double max_rt, Instrument::Type instrument_type, double resolution_ms1,
    double resolution_msn, double reference_mz, RawData::Polarity polarity);
}  // namespace XmlReader

#endif /* GRID_XMLREADER_HPP */
