#include "base64.hpp"

Base64::Base64(char *string_pointer, int precision, int little_endian)
    : m_string_pointer(string_pointer),
      m_precision(precision),
      m_little_endian(little_endian) {}

uint32_t Base64::get_uint32() {
    unsigned int b = 0;
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
    uint64_t b = 0;
    uint32_t b1 = get_uint32();
    uint32_t b2 = get_uint32();
    b = (uint64_t)b1 << 32;
    b |= b2;
    return b;
}

double Base64::get_double() {
    if (m_precision == 32) {
        uint32_t b = get_uint32();
        return *(float *)&b;
    } else if (m_precision == 64) {
        uint64_t b = get_uint64();
        return *(double *)&b;
    }
    return 0;
}
