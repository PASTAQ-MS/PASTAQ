#ifndef UTILS_COMPRESSION_HPP
#define UTILS_COMPRESSION_HPP

#include <zlib.h>
#include <iostream>
#include <streambuf>
#include <vector>

// This namespace contains necessary functions to (de)compress raw data.
namespace Compression {

enum state { OK, ERROR };

// Decompress raw data.
int inflate(std::vector<uint8_t> &in_data, std::vector<uint8_t> &out_data,
            size_t decompressed_len);

// Streambuf class allows a stream to write compressed data to a file by use of
// an intermediate buffer.
class DeflateStreambuf : public std::streambuf {
    // Buffer to store information before compression.
    char *buffer;
    size_t buffer_size;

    // File to write compressed data to.
    FILE *out_file = nullptr;

    // Zlib stream used in compression.
    z_stream strm;

   public:
    // Constructor sets buffer size.
    DeflateStreambuf(size_t _buffer_size = 16384);
    // Destructor flushes the buffer, closes the file, and deletes buffer.
    virtual ~DeflateStreambuf();

    // Open file and allocate Zlib state.
    int open(std::string const &filename);

   private:
    virtual int overflow(int c);  // Writes byte when buffer is full.
    virtual int sync();           // Flushes the buffer.
    int write_buffer(int flush);  // Compress data form buffer to file.
};

// DeflateStream uses the DeflateStreambuf to compress the data and write to a
// file.
class DeflateStream : private DeflateStreambuf, public std::ostream {
   public:
    DeflateStream(size_t buffer_size = 16384)
        : DeflateStreambuf(buffer_size), std::ostream(this) {}
    DeflateStream(std::string const &filename, size_t buffer_size = 16384)
        : DeflateStreambuf(buffer_size), std::ostream(this) {
        open(filename);
    }

    // Open streambuf and check for success.
    void open(std::string const &filename);
};

// Streambuf class allows a stream to read data from a file and decompress it
// using an intermediate buffer
class InflateStreambuf : public std::streambuf {
    // Buffer to store decompressed data.
    char *buffer;
    size_t buffer_size;

    // File to read compressed data from.
    FILE *in_file = nullptr;

    // Zlib stream used in decompression.
    z_stream strm;

   public:
    // Constructor sets buffer size.
    InflateStreambuf(size_t _buffer_size = 16384);
    // Destructor flushes the buffer, closes the file, and deletes buffer.
    virtual ~InflateStreambuf();

    // Open file and allocate Zlib state.
    int open(std::string const &filename);

   private:
    virtual int underflow();  // Read byte when buffer is empty.
    int read_buffer();        // Decompress data from file into the buffer.
};

// InflateStream uses the InflateStreambuf to decompress the data read from a
// file.
class InflateStream : private InflateStreambuf, public std::istream {
   public:
    InflateStream(size_t buffer_size = 16384)
        : InflateStreambuf(buffer_size), std::istream(this) {}
    InflateStream(std::string const &filename, size_t buffer_size = 16384)
        : InflateStreambuf(buffer_size), std::istream(this) {
        open(filename);
    }

    // Open streambuf and check for success.
    void open(std::string const &filename);
};

}  // namespace Compression

#endif /* UTILS_COMPRESSION_HPP */
