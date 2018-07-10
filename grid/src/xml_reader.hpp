#ifndef GRID_XMLREADER_HPP
#define GRID_XMLREADER_HPP

#include <map>
#include <optional>

#include "grid.hpp"

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

}  // namespace XmlReader

#endif /* GRID_XMLREADER_HPP */
