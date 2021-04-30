#include <zlib.h>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <vector>

#include "compression.hpp"

// decompress in chunks of 256 KB, this number can be increase or decreased
// based on available memory.
#define CHUNK 262144

// Decompress raw memory from the in_data vector into the out_data vector.
// Function allocates memory for the output vector. Function takes the length of
// the data after decompression.
int inflate(std::vector<uint8_t> &in_data, std::vector<uint8_t> &out_data,
            size_t decompressed_len) {
    int ret;
    z_stream strm;
    out_data.resize(decompressed_len);

    // allocate inflate state
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
    // decompress until deflate stream ends or end of file
    do {
        // read the amount of CHUNK or until the end of the data
        if (strm.total_in + CHUNK > in_data.size()) {
            strm.avail_in = in_data.size() - strm.total_in;
        } else {
            strm.avail_in = CHUNK;
        }
        if (strm.avail_in == 0) {  // end of data
            break;
        }
        // use the total bytes read as an offset for the input data
        strm.next_in = in + strm.total_in;

        // run inflate() on input until output buffer not full
        do {
            bytes_decompressed = strm.total_out;

            // calculate amount of free bytes in output buffer
            if (bytes_decompressed + CHUNK > decompressed_len) {
                strm.avail_out = decompressed_len - bytes_decompressed;
            } else {
                strm.avail_out = CHUNK;
            }

            // use next section of output buffer as output buffer
            strm.next_out = reinterpret_cast<unsigned char *>(
                &out_data[bytes_decompressed]);

            ret = inflate(&strm, Z_NO_FLUSH);
            assert(ret != Z_STREAM_ERROR);
            switch (ret) {
                case Z_NEED_DICT:
                    ret = Z_DATA_ERROR;
                case Z_DATA_ERROR:
                case Z_MEM_ERROR:
                    (void)inflateEnd(&strm);
                    return ret;
            }
        } while (bytes_decompressed == decompressed_len);

        // done when inflate() says it's done
    } while (ret != Z_STREAM_END);

    // clean up and return
    (void)inflateEnd(&strm);
    return ret == Z_STREAM_END ? Z_OK : Z_DATA_ERROR;
}
