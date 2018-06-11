#include <map>
#include <optional>

#include "grid.hpp"

namespace XmlReader {

struct Tag {
    std::string name;
    std::map<std::string, std::string> attributes;
    bool closed;
};

struct Scan {
    double mz;
    double rt;
    double value;
};

std::vector<Scan>
read_next_scan(std::istream& stream, Grid::Parameters& parameters);
std::optional<Tag> read_tag(std::istream& stream);

}  // namespace XmlReader
