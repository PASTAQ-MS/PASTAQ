// This namespace contains necessary functions to serialize commonly used types
// into a binary stream using the little endian byte order.
namespace Serialization {

// Write/read a single byte to/from the stream.
bool read_uint8(std::istream &stream, uint8_t *value);
bool write_uint8(std::ostream &stream, uint8_t value);

// Write/read an uint16 to/from the stream.
bool read_uint16(std::istream &stream, uint16_t *value);
bool write_uint16(std::ostream &stream, uint16_t value);

// Write/read an uint32 to/from the stream.
bool read_uint32(std::istream &stream, uint32_t *value);
bool write_uint32(std::ostream &stream, uint32_t value);

// Write/read an uint64 to/from the stream.
bool read_uint64(std::istream &stream, uint64_t *value);
bool write_uint64(std::ostream &stream, uint64_t value);

bool read_float(std::istream &stream, float *value);
bool write_float(std::ostream &stream, float value);

bool read_double(std::istream &stream, double *value);
bool write_double(std::ostream &stream, double value);

}  // namespace Serialization
