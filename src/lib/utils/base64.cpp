#include <cstring>

#include "base64.hpp"

inline uint32_t swap_uint32(const uint32_t &b, bool little_endian) {
    uint32_t ret;
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

Base64::Base64(unsigned char *string_pointer, int precision, bool little_endian)
    : m_string_pointer(string_pointer),
      m_precision(precision),
      m_little_endian(little_endian) {}

uint32_t Base64::get_uint32() {
    uint32_t b = 0;
    switch (m_bit) {
        case 0:
            b = m_translation_table[*m_string_pointer++] << 26;
            b |= m_translation_table[*m_string_pointer++] << 20;
            b |= m_translation_table[*m_string_pointer++] << 14;
            b |= m_translation_table[*m_string_pointer++] << 8;
            b |= m_translation_table[*m_string_pointer++] << 2;
            b |= m_translation_table[*m_string_pointer] >> 4;
            m_bit = 2;
            break;
        case 2:
            b = m_translation_table[*m_string_pointer++] << 28;
            b |= m_translation_table[*m_string_pointer++] << 22;
            b |= m_translation_table[*m_string_pointer++] << 16;
            b |= m_translation_table[*m_string_pointer++] << 10;
            b |= m_translation_table[*m_string_pointer++] << 4;
            b |= m_translation_table[*m_string_pointer] >> 2;
            m_bit = 4;
            break;
        case 4:
            b = m_translation_table[*m_string_pointer++] << 30;
            b |= m_translation_table[*m_string_pointer++] << 24;
            b |= m_translation_table[*m_string_pointer++] << 18;
            b |= m_translation_table[*m_string_pointer++] << 12;
            b |= m_translation_table[*m_string_pointer++] << 6;
            b |= m_translation_table[*m_string_pointer++];
            m_bit = 0;
            break;
    }
    return b;
}

uint64_t Base64::get_uint64() {
    uint64_t b1 = get_uint32();
    uint64_t b2 = get_uint32();
    b1 <<= 32;
    b1 |= b2;
    return b1;
}

double Base64::get_double() {
    if (m_precision == 32) {
        uint32_t b = get_uint32();
        b = swap_uint32(b, m_little_endian);
        float ret;
        std::memcpy(&ret, &b, sizeof(b));
        return ret;
    }
    if (m_precision == 64) {
        uint64_t b = get_uint64();
        b = swap_uint64(b, m_little_endian);
        double ret;
        std::memcpy(&ret, &b, sizeof(b));
        return ret;
    }
    return 0;
}
