#ifndef TESTUTILS_MOCKSTREAM_HPP
#define TESTUTILS_MOCKSTREAM_HPP

#include <iostream>
#include <streambuf>
#include <vector>

// Helper class to use a vector to mock istreams/ostreams. Note that the
// lifetime of the VectorStream class is not tied to the object we have mapped.
// If the vector get's deallocated or it's dimensions change, it will result in
// an error when using the stream. This is not a robust class and should be used
// with caution.
template <class T>
struct MockStream : public std::iostream {
    struct VectorStream : public std::streambuf {
        VectorStream(std::vector<T> &data) {
            auto begin = reinterpret_cast<char *>(&data[0]);
            auto end =
                reinterpret_cast<char *>(&data[0]) + sizeof(T) * data.size();
            setg(begin, begin, end);
            setp(begin, end);
        }
        pos_type seekoff(
            off_type off, std::ios_base::seekdir dir,
            std::ios_base::openmode which = std::ios_base::in) override {
            if (which == std::ios_base::in) {
                switch (dir) {
                    case std::ios::beg: {
                        setg(eback(), eback() + off, egptr());
                    } break;
                    case std::ios::cur: {
                        setg(eback(), gptr() + off, egptr());
                    } break;
                    case std::ios::end: {
                        setg(eback(), egptr() + off, egptr());
                    } break;
                    default: {
                        // ...
                    } break;
                }
            }
            return gptr() - eback();
        }
    } m_vs;
    MockStream(std::vector<T> &data) : m_vs(data), std::iostream(&m_vs) {}
    pos_type seekg(off_type off, std::ios_base::seekdir dir,
                   std::ios_base::openmode which = std::ios_base::in) {
        if (!good()) {
            return -1;
        }
        return m_vs.pubseekoff(off, dir, which);
    }
    pos_type tellg() {
        if (!good()) {
            return -1;
        }
        return m_vs.pubseekoff(0, std::ios::cur);
    }
};

#endif /* TESTUTILS_MOCKSTREAM_HPP */
