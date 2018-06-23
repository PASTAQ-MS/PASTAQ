namespace Endian {

constexpr bool platform_little_endian() {
    // NOTE(alex): This is a stub. Right now the code is designed to be used on
    // a x86_64 machine, which is LITTLE ENDIAN. With C++20 we will be able to
    // use std::endian::native == std::endian::little. Since currently no good
    // options exist to portably detect the machine endianness at compile time I
    // will defer the implementation of this function to the future.
    return true;
}

inline uint32_t swap_uint32(uint32_t b, bool little_endian) {
    if (!platform_little_endian() && !little_endian) {
        return (b >> 24 & 0xFFFFFFFF) | (b >> 8 & 0x0000FF00) |
               (b << 8 & 0x00FF0000) | (b << 24 & 0xFF000000);
    }
    if (platform_little_endian() && little_endian) {
        return (b << 24 & 0xFFFFFFFF) | (b << 8 & 0x0000FF00) |
               (b >> 8 & 0x00FF0000) | (b >> 24 & 0xFF000000);
    }
    return b;
}

inline uint64_t swap_uint64(uint64_t b, bool little_endian) {
    if (!platform_little_endian() && !little_endian) {
        uint64_t b1 = (b >> 32) & 0x00000000FFFFFFFF;
        uint64_t b2 = b & 0x00000000FFFFFFFF;
        b1 = swap_uint32(b1, little_endian);
        b2 = swap_uint32(b2, little_endian);
        return (b2 << 32) | b1;
    }
    if (platform_little_endian() && little_endian) {
        uint64_t b1 = b & 0x00000000FFFFFFFF;
        uint64_t b2 = (b >> 32) & 0x00000000FFFFFFFF;
        b1 = swap_uint32(b1, little_endian);
        b2 = swap_uint32(b2, little_endian);
        return (b1 << 32) | b2;
    }
    return b;
}

}  // namespace Endian
