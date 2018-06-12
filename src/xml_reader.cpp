#include <iostream>
#include <regex>
#include <sstream>

#include "xml_reader.hpp"

// Read the next mz1 scan from the stream.
std::optional<std::vector<XmlReader::Scan>> XmlReader::read_next_scan(
    std::istream& stream, Grid::Parameters& parameters) {
    while (stream.good()) {
        auto tag = XmlReader::read_tag(stream);
        if (!tag) {
            continue;
        }
        if (tag.value().name == "scan" && !tag.value().closed) {
            auto scan_attributes = tag.value().attributes;
            // We are only interested on ms1 scans for now.
            if (scan_attributes["msLevel"] != "1") {
                continue;
            }
            // Extract the peak count.
            if (scan_attributes.find("peaksCount") == scan_attributes.end()) {
                return std::nullopt;
            }
            int peaks_count = std::stoi(scan_attributes["peaksCount"]);

            // Extract the retention time.
            if (scan_attributes.find("retentionTime") ==
                scan_attributes.end()) {
                return std::nullopt;
            }
            std::regex rt_regex("PT([[:digit:]]+.?[[:digit:]]+)S");
            std::smatch matches;
            if (!std::regex_search(scan_attributes["retentionTime"], matches,
                                   rt_regex) ||
                matches.size() != 2) {
                return std::nullopt;
            }
            double retention_time = std::stod(matches[1]);

            // Check if we are on the desired region as defined by Grid::Bounds.
            if (retention_time < parameters.bounds.min_rt) {
                continue;
            }
            // Assuming linearity of the retention time on the mzXML file. We
            // stop searching for the scan, since we are out of bounds.
            if (retention_time > parameters.bounds.max_rt) {
                return std::nullopt;
            }

            // Fetch the peaks tag.
            tag = XmlReader::read_tag(stream);
            while (tag) {
                if (tag.value().name == "scan" && tag.value().closed) {
                    return std::nullopt;
                }
                if (tag.value().name == "peaks") {
                    break;
                }
                tag = XmlReader::read_tag(stream);
            }
            auto peak_attributes = tag.value().attributes;

            // Extract the precision from the peaks tag.
            if (peak_attributes.find("precision") == peak_attributes.end()) {
                return std::nullopt;
            }
            int precision = std::stoi(peak_attributes["precision"]);

            // Extract the byteOrder from the peaks tag. This determines the
            // endianness in which the data was stored. `network` ==
            // `big_endian`.
            if (peak_attributes.find("byteOrder") == peak_attributes.end()) {
                return std::nullopt;
            }
            auto byte_order = peak_attributes["byteOrder"];
            auto little_endian = (byte_order == "network") ? false : true;
            std::cout << little_endian << std::endl;

            // Extract the contentType/pairOrder from the peaks tag and exit if
            // it is not `m/z-int`. In older versions of ProteoWizard, the
            // conversion was not validated and the tag was incorrect. Here we
            // are supporting both versions for compatibility but we are not
            // trying to be exhaustive.
            auto content_type_found =
                peak_attributes.find("contentType") != peak_attributes.end();
            auto pair_order_found =
                peak_attributes.find("pairOrder") != peak_attributes.end();
            if (!content_type_found && !pair_order_found) {
                return std::nullopt;
            }
            std::string pair_order;
            if (content_type_found) {
                pair_order = peak_attributes["contentType"];
            } else {
                pair_order = peak_attributes["pairOrder"];
            }

            // TODO: Extract the peaks from the data, performing the Base64
            // decoding in the process.
            // TODO: We don't need to extract the peaks when we are not inside
            // the mz bounds.
 
            std::cout << "OKI" << std::endl;
            return std::nullopt;
            std::vector<XmlReader::Scan> scans;
            return scans;
        }
    }
    return std::nullopt;
}

// TODO: Should fix the parser to accept spaces in quoted strings.
// TODO: This is a very naive way of performing the tag parsing, going through
// it character by character. We should evaluate the performance and see if it
// is worth to make it faster by buffering.
std::optional<XmlReader::Tag> XmlReader::read_tag(std::istream& stream) {
    bool is_closed = false;
    bool reading_content = false;
    auto tag = std::optional<Tag>(Tag{});

    // Store the tag contents in a buffer for further processing.
    std::string buffer;
    while (stream.good()) {
        char c = stream.get();
        if (c == '<') {
            if (stream.peek() == '/') {
                tag->closed = true;
                stream.get();
            }
            reading_content = true;
        } else if (c == '>') {
            is_closed = true;
            break;
        } else if (reading_content) {
            buffer += c;
        }
    }

    if (!is_closed) {
        return std::nullopt;
    }

    // Tokenize tag contents.
    std::stringstream ss(buffer);
    std::string item;
    std::vector<std::string> tokens;
    while (ss >> item) {
        tokens.push_back(item);
    };

    // Fill tag fieds.
    for (int i = 0; i < tokens.size(); ++i) {
        auto token = tokens[i];
        if (i == 0) {
            tag->name = token;
        } else {
            // Parse each attribute on the hash map.
            std::regex attribute_regex("(\\S+)=\"(\\S+)\"");
            std::smatch matches;
            if (std::regex_search(token, matches, attribute_regex)) {
                // Malformed attributes, expected `name=value`
                if (matches.size() != 3 || matches[1] == "" ||
                    matches[2] == "") {
                    return std::nullopt;
                }
                tag->attributes[matches[1]] = matches[2];
            } else {
                return std::nullopt;
            }
        }
    }
    return tag;
};
