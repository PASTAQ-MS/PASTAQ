#ifndef GRID_XMLREADER_HPP
#define GRID_XMLREADER_HPP

#include <map>
#include <optional>

#include "grid.hpp"

// TODO: Add documentation for this namespace
namespace XmlReader {

struct Tag {
    std::string name;
    std::map<std::string, std::string> attributes;
    bool closed;
};

std::optional<std::vector<Grid::Peak>> read_next_scan(
    std::istream &stream, Grid::Parameters &parameters);
std::optional<Tag> read_tag(std::istream &stream);

// Read data until the next tag is found and trim whitespace at the beginning in
// necessary.
std::optional<std::string> read_data(std::istream &stream);

}  // namespace XmlReader

namespace Base64 {

const unsigned char translation_table[256] = {
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 62, 64, 64, 64, 63, 52, 53, 54, 55, 56, 57, 58, 59, 60,
    61, 64, 64, 64, 64, 64, 64, 64, 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10,
    11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 64, 64, 64, 64,
    64, 64, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42,
    43, 44, 45, 46, 47, 48, 49, 50, 51, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64};

inline static uint32_t get32(char *&p, int &bit, bool little_endian) {
    // bit is either 0, 2, 4
    // indicates how many bits should be discarded
    unsigned int b;

    switch (bit) {
        case 0:
            b = translation_table[*p++] << 26;
            b |= translation_table[*p++] << 20;
            b |= translation_table[*p++] << 14;
            b |= translation_table[*p++] << 8;
            b |= translation_table[*p++] << 2;
            b |= translation_table[*p] >> 4;
            bit = 2;
            break;
        case 2:
            b = translation_table[*p++] << 28;
            b |= translation_table[*p++] << 22;
            b |= translation_table[*p++] << 16;
            b |= translation_table[*p++] << 10;
            b |= translation_table[*p++] << 4;
            b |= translation_table[*p] >> 2;
            bit = 4;
            break;
        case 4:
            b = translation_table[*p++] << 30;
            b |= translation_table[*p++] << 24;
            b |= translation_table[*p++] << 18;
            b |= translation_table[*p++] << 12;
            b |= translation_table[*p++] << 6;
            b |= translation_table[*p++];
            bit = 0;
            break;
    }
    // Swap::MakeInt32(b, little_endian);
    return b;
}

inline static uint64_t get64(char *&p, int &bit, bool little_endian) {
    // bit is either 0, 2, 4
    // indicates how many bits should be discarded
    uint64_t b;
    uint32_t b1, b2;

    switch (bit) {
        case 0:
            b1 = translation_table[*p++] << 26;
            b1 |= translation_table[*p++] << 20;
            b1 |= translation_table[*p++] << 14;
            b1 |= translation_table[*p++] << 8;
            b1 |= translation_table[*p++] << 2;
            b1 |= translation_table[*p] >> 4;
            bit = 2;
            break;
        case 2:
            b1 = translation_table[*p++] << 28;
            b1 |= translation_table[*p++] << 22;
            b1 |= translation_table[*p++] << 16;
            b1 |= translation_table[*p++] << 10;
            b1 |= translation_table[*p++] << 4;
            b1 |= translation_table[*p] >> 2;
            bit = 4;
            break;
        case 4:
            b1 = translation_table[*p++] << 30;
            b1 |= translation_table[*p++] << 24;
            b1 |= translation_table[*p++] << 18;
            b1 |= translation_table[*p++] << 12;
            b1 |= translation_table[*p++] << 6;
            b1 |= translation_table[*p++];
            bit = 0;
            break;
    }

    switch (bit) {
        case 0:
            b2 = translation_table[*p++] << 26;
            b2 |= translation_table[*p++] << 20;
            b2 |= translation_table[*p++] << 14;
            b2 |= translation_table[*p++] << 8;
            b2 |= translation_table[*p++] << 2;
            b2 |= translation_table[*p] >> 4;
            bit = 2;
            break;
        case 2:
            b2 = translation_table[*p++] << 28;
            b2 |= translation_table[*p++] << 22;
            b2 |= translation_table[*p++] << 16;
            b2 |= translation_table[*p++] << 10;
            b2 |= translation_table[*p++] << 4;
            b2 |= translation_table[*p] >> 2;
            bit = 4;
            break;
        case 4:
            b2 = translation_table[*p++] << 30;
            b2 |= translation_table[*p++] << 24;
            b2 |= translation_table[*p++] << 18;
            b2 |= translation_table[*p++] << 12;
            b2 |= translation_table[*p++] << 6;
            b2 |= translation_table[*p++];
            bit = 0;
            break;
    }
    b = (uint64_t)b1 << 32;
    b |= b2;
    // Swap::MakeInt64(b, little_endian);
    return b;
}

// TODO: Can we make it so that we don't have to pass the 'bit' variable around?
inline static std::optional<double> get_double(char *&p, int &bit,
                                               int precision,
                                               bool little_endian) {
    if (precision == 32) {
        uint32_t b = get32(p, bit, little_endian);
        return *(float *)&b;
    } else if (precision == 64) {
        uint64_t b = get64(p, bit, little_endian);
        return *(double *)&b;
    }
    return std::nullopt;
}

}  // namespace Base64

#endif /* GRID_XMLREADER_HPP */
