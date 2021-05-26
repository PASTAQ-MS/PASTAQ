#include <zlib.h>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <vector>

#include "compression.hpp"

// Decompress in chunks of 256 KB, this number can be increase or decreased
// based on available memory.
#define CHUNK 262144

// Decompress raw memory from the in_data vector into the out_data vector.
// Function allocates memory for the output vector. Function takes the length of
// the data after decompression.
int Compression::inflate(std::vector<uint8_t> &in_data,
                         std::vector<uint8_t> &out_data,
                         size_t decompressed_len) {
    int ret;
    z_stream strm;
    out_data.resize(decompressed_len);

    // Allocate inflate state.
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    strm.avail_in = 0;
    strm.next_in = Z_NULL;
    ret = inflateInit(&strm);
    if (ret != Z_OK) {
        return ret;
    }

    unsigned char *in = reinterpret_cast<unsigned char *>(&in_data[0]);
    size_t bytes_decompressed = 0;
    // Decompress until deflate stream ends or end of file.
    do {
        // Read the amount of CHUNK or until the end of the data.
        if (strm.total_in + CHUNK > in_data.size()) {
            strm.avail_in = in_data.size() - strm.total_in;
        } else {
            strm.avail_in = CHUNK;
        }
        if (strm.avail_in == 0) {  // End of data.
            break;
        }
        // Use the total bytes read as an offset for the input data.
        strm.next_in = in + strm.total_in;

        // Run inflate() on input until output buffer not full.
        do {
            bytes_decompressed = strm.total_out;

            // Calculate amount of free bytes in output buffer.
            if (bytes_decompressed + CHUNK > decompressed_len) {
                strm.avail_out = decompressed_len - bytes_decompressed;
            } else {
                strm.avail_out = CHUNK;
            }

            // Use next section of output buffer as output buffer.
            strm.next_out = reinterpret_cast<unsigned char *>(
                &out_data[bytes_decompressed]);

            ret = inflate(&strm, Z_NO_FLUSH);
            assert(ret != Z_STREAM_ERROR);

            switch (ret) {
                case Z_NEED_DICT:
                    ret = Z_DATA_ERROR;
                    (void)inflateEnd(&strm);
                    return ret;
                case Z_DATA_ERROR:
                    (void)inflateEnd(&strm);
                    return ret;
                case Z_MEM_ERROR:
                    (void)inflateEnd(&strm);
                    return ret;
            }
        } while (bytes_decompressed == decompressed_len);

        // Done when inflate() says it's done.
    } while (ret != Z_STREAM_END);

    // Clean up and return.
    (void)inflateEnd(&strm);
    return ret == Z_STREAM_END ? Z_OK : Z_DATA_ERROR;
}

// Initialize buffer.
Compression::DeflateStreambuf::DeflateStreambuf(size_t _buffer_size)
    : buffer_size(_buffer_size) {
    // Allocate new buffer.
    buffer = new char[buffer_size];
    // Set streambuf's internal pointer to buffer.
    setp(buffer, buffer + buffer_size);
}

// Destructor flushes the buffer, deletes the allocated memory, closes the
// output file, and frees the allocated Zlib state.
Compression::DeflateStreambuf::~DeflateStreambuf() {
    // Flush current buffer.
    sync();

    // Perform last write for Zlib to flush it's internal buffer.
    strm.next_in = nullptr;
    strm.avail_in = 0;
    write_buffer(Z_FINISH);

    delete[] buffer;
    if (out_file) {
        fclose(out_file);
    }
    (void)deflateEnd(&strm);
}

// Open file, allocate buffer, and initialize Zlib state.
int Compression::DeflateStreambuf::open(std::string const &filename) {
    // Open file.
    out_file = fopen(filename.c_str(), "wb");
    if (out_file == NULL) {
        return ERROR;
    }

    // Initialize zlib stream.
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    int ret = deflateInit(&strm, Z_DEFAULT_COMPRESSION);
    if (ret != Z_OK) {
        return ERROR;
    }
    return OK;
}

// This function writes a character when the buffer is full.
int Compression::DeflateStreambuf::overflow(int c) {
    // Flush buffer.
    if (sync() == -1) {
        // Signal error when sync fails
        return EOF;
    }

    // Set pointer to buffer.
    setp(buffer, buffer + buffer_size);
    return sputc(c);
}

// Flushes the buffer using the Zlib library for compression.
int Compression::DeflateStreambuf::sync() {
    if (pptr() > pbase()) {  // buffer not empty
        int ret = write_buffer(Z_NO_FLUSH);
        if (ret != Z_OK) {
            return -1;
        }
        // Set pointer to buffer.
        setp(buffer, buffer + buffer_size);
    }
    return 0;
}

// Write the current buffer using the Zlib library. While this flushes our
// intermediate buffer, the Zlib library may not write all of the compressed
// data to the file immediately.
// When writing the buffer in between compression, the argument flush should
// be equal to Z_NO_FINISH. After writing the final block of data, the
// argument flush should be equal to Z_FINISH so Zlib knows to flush all of
// the compressed data.
int Compression::DeflateStreambuf::write_buffer(int flush) {
    // Set buffer and buffer size.
    strm.avail_in = pptr() - pbase();
    strm.next_in = reinterpret_cast<unsigned char *>(pbase());

    int ret;
    // Create and set intermediate buffer for compressed data.
    unsigned char *out = new unsigned char[buffer_size];
    // Deflate and write compressed data to file until we do not fill the
    // output buffer.
    do {
        strm.avail_out = buffer_size;
        strm.next_out = out;

        // Deflate buffer.
        ret = deflate(&strm, flush);

        // Write compressed data to file.
        size_t have = buffer_size - strm.avail_out;

        if (fwrite(out, 1, have, out_file) != have || ferror(out_file)) {
            return Z_ERRNO;
        }
    } while (strm.avail_out == 0);

    delete[] out;
    return ret;
}

// Open streambuf and check for success.
void Compression::DeflateStream::open(std::string const &filename) {
    int state = DeflateStreambuf::open(filename);
    if (state == ERROR) {
        setstate(std::ios::badbit);
    }
}

// Initialize buffer.
Compression::InflateStreambuf::InflateStreambuf(size_t _buffer_size)
    : buffer_size(_buffer_size) {
    // Allocate new buffer.
    buffer = new char[buffer_size];
    // Set streambuf's internal pointer to buffer.
    setg(buffer, buffer + buffer_size, buffer + buffer_size);
}

// Destructor deleted the allocated memory, closes the input file, and frees the
// allocated Zlib state.
Compression::InflateStreambuf::~InflateStreambuf() {
    delete[] buffer;
    if (in_file) {
        fclose(in_file);
    }
    (void)inflateEnd(&strm);
}

// Open file and initialize Zlib state.
int Compression::InflateStreambuf::open(std::string const &filename) {
    // Open file.
    in_file = fopen(filename.c_str(), "rb");
    if (in_file == NULL) {
        return ERROR;
    }

    // Initialize zlib stream.
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;
    int ret = inflateInit(&strm);
    if (ret != Z_OK) {
        return ERROR;
    }
    return OK;
}

// Read a character when the buffer is empty.
int Compression::InflateStreambuf::underflow() {
    if (gptr() < egptr()) {
        return *gptr();
    }

    size_t nread = read_buffer();

    if (nread <= 0) {
        return EOF;
    }
    setg(buffer, buffer, buffer + nread);
    return static_cast<unsigned char>(*gptr());
}

// Fill the buffer using the Zlib library for decompression, returns the number
// of bytes that were read to the buffer.
int Compression::InflateStreambuf::read_buffer() {
    // Buffer to store data read from file.
    unsigned char *in = new unsigned char[buffer_size];
    strm.avail_in = fread(in, 1, buffer_size, in_file);
    if (ferror(in_file)) {
        (void)inflateEnd(&strm);
        return 0;
    }

    // Return on EOF.
    if (strm.avail_in == 0) {
        return 0;
    }
    strm.next_in = in;

    size_t bytes_read = 0;
    // Read and inflate data until the output buffer is full or the input buffer
    // is empty.
    do {
        // Set output buffer.
        size_t size = buffer_size - bytes_read;
        strm.avail_out = size;
        strm.next_out = reinterpret_cast<unsigned char *>(eback() + bytes_read);

        // Inflate data.
        int ret = inflate(&strm, Z_NO_FLUSH);
        switch (ret) {
            case Z_NEED_DICT:
                ret = Z_DATA_ERROR;
                [[fallthrough]];
            case Z_DATA_ERROR:
            case Z_MEM_ERROR:
                (void)inflateEnd(&strm);
                return ret;
        }
        bytes_read += size - strm.avail_out;
    } while (strm.avail_in != 0 && strm.avail_out != 0);

    // Move to position in file directly after the last byte that was
    // decompressed already.
    fseek(in_file, -(int)strm.avail_in, SEEK_CUR);

    delete[] in;
    return bytes_read;
}

// Open streambuf and check for success.
void Compression::InflateStream::open(std::string const &filename) {
    int state = InflateStreambuf::open(filename);
    if (state == ERROR) {
        setstate(std::ios::badbit);
    }
}
