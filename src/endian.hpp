#ifndef GRID_ENDIAN_HPP
#define GRID_ENDIAN_HPP

namespace Endian {

inline uint32_t swap_uint32(const uint32_t &b, bool little_endian) {
    uint64_t ret;
    uint8_t *data = reinterpret_cast<uint8_t *>(&ret);

    if (little_endian) {
        data[0] = b >> 24;
        data[1] = b >> 16;
        data[2] = b >> 8;
        data[3] = b >> 0;
    } else {
        data[3] = b >> 24;
        data[2] = b >> 16;
        data[1] = b >> 8;
        data[0] = b >> 0;
    }
    return ret;
}

inline uint64_t swap_uint64(const uint64_t &b, bool little_endian) {
    uint64_t ret;
    uint8_t *data = reinterpret_cast<uint8_t *>(&ret);

    if (little_endian) {
        data[0] = b >> 56;
        data[1] = b >> 48;
        data[2] = b >> 40;
        data[3] = b >> 32;
        data[4] = b >> 24;
        data[5] = b >> 16;
        data[6] = b >> 8;
        data[7] = b >> 0;
    } else {
        data[7] = b >> 56;
        data[6] = b >> 48;
        data[5] = b >> 40;
        data[4] = b >> 32;
        data[3] = b >> 24;
        data[2] = b >> 16;
        data[1] = b >> 8;
        data[0] = b >> 0;
    }
    return ret;
}

}  // namespace Endian

#endif /* GRID_ENDIAN_HPP */
