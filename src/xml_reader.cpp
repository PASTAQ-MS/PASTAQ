#include <iostream>
#include <regex>
#include <sstream>

#include "xml_reader.hpp"

// Read the next mz1 scan from the stream.
std::vector<XmlReader::Scan> XmlReader::read_next_scan(
    std::istream& stream, Grid::Parameters& parameters) {
    std::vector<XmlReader::Scan> scans;
    while (stream.good()) {
        auto tag = XmlReader::read_tag(stream);
        if (!tag) {
            continue;
        }
        if (tag.value().name == "scan" && !tag.value().closed) {
            std::cout << "SCAN FOUND:" << std::endl;
            for (const auto& e : tag.value().attributes) {
                std::cout << e.first << " = " << e.second << std::endl;
            }
            
        }
    }
    return scans;
}

// TODO: Should fix the parser to accept spaces in quoted strings.
// TODO: This is a very naive way of performing the tag parsing, going through it
// character by character. We should evaluate the performance and see if it is
// worth to make it faster by buffering.
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
