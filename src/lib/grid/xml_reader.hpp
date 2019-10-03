#ifndef GRID_XMLREADER_HPP
#define GRID_XMLREADER_HPP

#include <map>
#include <optional>

#include "grid/grid.hpp"
#include "grid/raw_data.hpp"

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
    double resolution_msn, double reference_mz, RawData::Polarity polarity,
    size_t ms_level);

// Read an entire mzIdentML file into a IdentData::IdentData data structure.
IdentData::IdentData read_mzidentml(std::istream &stream);
}  // namespace XmlReader

#endif /* GRID_XMLREADER_HPP */
