#include <cstring>
#include <string>
#include <vector>

#include "utils/base64.hpp"

void Base64::decode_base64(const std::string &input, std::vector<uint8_t> &output) {
    size_t in_len = input.size();

    size_t out_len = in_len / 4 * 3;
    if (input[in_len - 1] == '=') {
        --out_len;
    }
    if (input[in_len - 2] == '=') {
        --out_len;
    }

    output.resize(out_len);

    for (size_t i = 0, j = 0; i < in_len; i += 4, j += 3) {
        uint8_t a = input[i] == '='
                        ? 0
                        : translation_table[static_cast<int>(input[i])];
        uint8_t b = input[i + 1] == '='
                        ? 0
                        : translation_table[static_cast<int>(input[i + 1])];
        uint8_t c = input[i + 2] == '='
                        ? 0
                        : translation_table[static_cast<int>(input[i + 2])];
        uint8_t d = input[i + 3] == '='
                        ? 0
                        : translation_table[static_cast<int>(input[i + 3])];

        if (j < out_len) {
            output[j] = (a << 2) + (b >> 4);
        }
        if (j + 1 < out_len) {
            output[j + 1] = (b << 4) + (c >> 2);
        }
        if (j + 2 < out_len) {
            output[j + 2] = (c << 6) + d;
        }
    }
}

// Interpreting functions.

// Read four bytes from the stream from the start index, order bytes
// based on byte order.
uint32_t Base64::interpret_uint32(std::vector<uint8_t> &data, size_t offset,
                          bool little_endian) {
    if (data.size() < offset + 4) {
        return uint32_t{};
    }

    uint32_t ret;
    uint8_t *bytes = reinterpret_cast<uint8_t *>(&ret);
    if (little_endian) {
        bytes[0] = data[offset];
        bytes[1] = data[offset + 1];
        bytes[2] = data[offset + 2];
        bytes[3] = data[offset + 3];
    } else {
        bytes[0] = data[offset + 3];
        bytes[1] = data[offset + 2];
        bytes[2] = data[offset + 1];
        bytes[3] = data[offset];
    }
    return ret;
}

uint64_t Base64::interpret_uint64(std::vector<uint8_t> &data, size_t offset,
                          bool little_endian) {
    if (data.size() < offset + 8) {
        return uint64_t{};
    }

    uint64_t ret;
    uint8_t *bytes = reinterpret_cast<uint8_t *>(&ret);
    if (little_endian) {
        bytes[0] = data[offset];
        bytes[1] = data[offset + 1];
        bytes[2] = data[offset + 2];
        bytes[3] = data[offset + 3];
        bytes[4] = data[offset + 4];
        bytes[5] = data[offset + 5];
        bytes[6] = data[offset + 6];
        bytes[7] = data[offset + 7];
    } else {
        bytes[0] = data[offset + 7];
        bytes[1] = data[offset + 6];
        bytes[2] = data[offset + 5];
        bytes[3] = data[offset + 4];
        bytes[4] = data[offset + 3];
        bytes[5] = data[offset + 2];
        bytes[6] = data[offset + 1];
        bytes[7] = data[offset];
    }
    return ret;
}

// Returns the float represented by the data vector at offset, interpreted using
// the specified byte order.
float Base64::interpret_float(std::vector<uint8_t> &data, size_t offset,
                      bool little_endian) {
    uint32_t bytes = interpret_uint32(data, offset, little_endian);
    float ret;
    std::memcpy(&ret, &bytes, sizeof(bytes));
    return ret;
}

// Returns the double represented by the data vector at offset, interpreted
// using the specified byte order.
double Base64::interpret_double(std::vector<uint8_t> &data, size_t offset,
                        bool little_endian) {
    uint64_t bytes = interpret_uint64(data, offset, little_endian);
    double ret;
    std::memcpy(&ret, &bytes, sizeof(bytes));
    return ret;
}
