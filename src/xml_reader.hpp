#include <map>
#include <optional>

#include "grid.hpp"

namespace XmlReader {

struct Tag {
    std::string name;
    std::map<std::string, std::string> attributes;
    bool closed;
};

struct Scan {
    double mz;
    double rt;
    double value;
};

std::optional<std::vector<Scan>> read_next_scan(std::istream& stream,
                                                Grid::Parameters& parameters);
std::optional<Tag> read_tag(std::istream& stream);

// Read data until the next tag is found and trim whitespace at the beginning in
// necessary.
std::optional<std::string> read_data(std::istream& stream);

}  // namespace XmlReader

namespace Base64 {

/*
    inline static void getFloatFloat(char *&p, int &bit, double &f, double &i,
                                     int precision, int littleendian) {
        if (precision == 32) {
            LCMSInt32 b;
            b = get32(p, bit, littleendian);
            f = *(float *)&b;
            b = get32(p, bit, littleendian);
            i = *(float *)&b;
        } else if (precision == 64) {
            LCMSInt64 b;
            b = get64(p, bit, littleendian);
            f = *(double *)&b;
            b = get64(p, bit, littleendian);
            i = *(double *)&b;
        } else {
            std::cout << "Not handled precision!";
            exit(2);
        }
    }

    ...

    inline static unsigned LCMSInt32 get32(char *&p, int &bit,
                                           int littleendian) {
        // bit is either 0, 2, 4
        // indicates how many bits should be discarded
        unsigned int b;

        switch (bit) {
            case 0:
                b = pr2six[*p++] << 26;
                b |= pr2six[*p++] << 20;
                b |= pr2six[*p++] << 14;
                b |= pr2six[*p++] << 8;
                b |= pr2six[*p++] << 2;
                b |= pr2six[*p] >> 4;
                bit = 2;
                break;
            case 2:
                b = pr2six[*p++] << 28;
                b |= pr2six[*p++] << 22;
                b |= pr2six[*p++] << 16;
                b |= pr2six[*p++] << 10;
                b |= pr2six[*p++] << 4;
                b |= pr2six[*p] >> 2;
                bit = 4;
                break;
            case 4:
                b = pr2six[*p++] << 30;
                b |= pr2six[*p++] << 24;
                b |= pr2six[*p++] << 18;
                b |= pr2six[*p++] << 12;
                b |= pr2six[*p++] << 6;
                b |= pr2six[*p++];
                bit = 0;
                break;
        }
        Swap::MakeInt32(b, littleendian);
        return b;
    }

    inline static unsigned LCMSInt64 get64(char *&p, int &bit,
                                           int littleendian) {
        // bit is either 0, 2, 4
        // indicates how many bits should be discarded
        unsigned LCMSInt64 b;
        unsigned LCMSInt32 b1, b2;

        switch (bit) {
            case 0:
                b1 = pr2six[*p++] << 26;
                b1 |= pr2six[*p++] << 20;
                b1 |= pr2six[*p++] << 14;
                b1 |= pr2six[*p++] << 8;
                b1 |= pr2six[*p++] << 2;
                b1 |= pr2six[*p] >> 4;
                bit = 2;
                break;
            case 2:
                b1 = pr2six[*p++] << 28;
                b1 |= pr2six[*p++] << 22;
                b1 |= pr2six[*p++] << 16;
                b1 |= pr2six[*p++] << 10;
                b1 |= pr2six[*p++] << 4;
                b1 |= pr2six[*p] >> 2;
                bit = 4;
                break;
            case 4:
                b1 = pr2six[*p++] << 30;
                b1 |= pr2six[*p++] << 24;
                b1 |= pr2six[*p++] << 18;
                b1 |= pr2six[*p++] << 12;
                b1 |= pr2six[*p++] << 6;
                b1 |= pr2six[*p++];
                bit = 0;
                break;
        }

        switch (bit) {
            case 0:
                b2 = pr2six[*p++] << 26;
                b2 |= pr2six[*p++] << 20;
                b2 |= pr2six[*p++] << 14;
                b2 |= pr2six[*p++] << 8;
                b2 |= pr2six[*p++] << 2;
                b2 |= pr2six[*p] >> 4;
                bit = 2;
                break;
            case 2:
                b2 = pr2six[*p++] << 28;
                b2 |= pr2six[*p++] << 22;
                b2 |= pr2six[*p++] << 16;
                b2 |= pr2six[*p++] << 10;
                b2 |= pr2six[*p++] << 4;
                b2 |= pr2six[*p] >> 2;
                bit = 4;
                break;
            case 4:
                b2 = pr2six[*p++] << 30;
                b2 |= pr2six[*p++] << 24;
                b2 |= pr2six[*p++] << 18;
                b2 |= pr2six[*p++] << 12;
                b2 |= pr2six[*p++] << 6;
                b2 |= pr2six[*p++];
                bit = 0;
                break;
        }
        b = (unsigned LCMSInt64)b1 << 32;
        b |= b2;
        Swap::MakeInt64(b, littleendian);
        return b;
    }

    ...

    #include "Mesh/Base64.h"
    const unsigned char Base64::pr2six[256] =
    {
        64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
        64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
        64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 62, 64, 64, 64, 63,
        52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 64, 64, 64, 64, 64, 64,
        64,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
        15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 64, 64, 64, 64, 64,
        64, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
        41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 64, 64, 64, 64, 64,
        64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
        64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
        64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
        64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
        64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
        64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
        64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
        64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64
    };
*/

}  // namespace Base64
