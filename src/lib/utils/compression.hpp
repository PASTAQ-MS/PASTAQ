#ifndef UTILS_COMPRESSION_HPP
#define UTILS_COMPRESSION_HPP

#include <vector>

// This namespace contains necessary functions to (de)compress raw data.
namespace Compression {

// (De)compression.
int inflate(std::vector<uint8_t> &in_data, std::vector<uint8_t> &out_data,
            size_t decompressed_len);

}  // namespace Compression

#endif /* UTILS_COMPRESSION_HPP */
