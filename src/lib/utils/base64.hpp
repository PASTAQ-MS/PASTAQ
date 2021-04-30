#ifndef UTILS_BASE64_HPP
#define UTILS_BASE64_HPP

#include <vector>

// This namespace contains functions to decode base64-encoded data into raw
// memory.
namespace Base64 {

// Translation table used in decoding function.
static constexpr unsigned char translation_table[] = {
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

// Decode input string containing base64-encoded data into bytes, result is
// returned in the output vector. This function allocates the appropriate amount
// of memory for the output vector.
void decode_base64(std::string const &input, std::vector<uint8_t> &output);

// Interpreting raw data into floating-point values.
uint32_t interpret_uint32(std::vector<uint8_t> &data, size_t offset,
                          bool little_endian);
uint64_t interpret_uint64(std::vector<uint8_t> &data, size_t offset,
                          bool little_endian);
float interpret_float(std::vector<uint8_t> &data, size_t offset,
                      bool little_endian);
double interpret_double(std::vector<uint8_t> &data, size_t offset,
                        bool little_endian);

}  // namespace Base64

#endif /* UTILS_BASE64_HPP */
