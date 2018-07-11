#include "doctest.h"

#include "testutils/mock_stream.hpp"
#include "utils/serialization.hpp"

TEST_CASE("Read/Write double values") {
    std::vector<double> source_data = {1.0, 2.0, 3.0, 4.0, -1.0, -2.0, -3.0};
    std::vector<double> destination_data;
    std::vector<char> stream_data(sizeof(double) * source_data.size());
    MockStream<char> stream(stream_data);
    // Write data to stream.
    for (const auto &value : source_data) {
        CHECK(Serialization::write_double(stream, value));
    }
    // Read data from stream.
    for (const auto &value : source_data) {
        double read_value;
        CHECK(Serialization::read_double(stream, &read_value));
        CHECK(read_value == value);
    }
}

TEST_CASE("Read/Write float values") {
    std::vector<float> source_data = {1.0f,  2.0f,  3.0f, 4.0f,
                                      -1.0f, -2.0f, -3.0f};
    std::vector<float> destination_data;
    std::vector<char> stream_data(sizeof(float) * source_data.size());
    MockStream<char> stream(stream_data);
    // Write data to stream.
    for (const auto &value : source_data) {
        CHECK(Serialization::write_float(stream, value));
    }
    // Read data from stream.
    for (const auto &value : source_data) {
        float read_value;
        CHECK(Serialization::read_float(stream, &read_value));
        CHECK(read_value == value);
    }
}

TEST_CASE("Read/Write uint64_t values") {
    std::vector<uint64_t> source_data = {
        0x8877665544332211, 0x1122334455667788, 0, 1, 2, 3, 4};
    std::vector<uint64_t> destination_data;
    std::vector<char> stream_data(sizeof(uint64_t) * source_data.size());
    MockStream<char> stream(stream_data);
    // Write data to stream.
    for (const auto &value : source_data) {
        CHECK(Serialization::write_uint64(stream, value));
    }
    // Read data from stream.
    for (const auto &value : source_data) {
        uint64_t read_value;
        CHECK(Serialization::read_uint64(stream, &read_value));
        CHECK(read_value == value);
    }
}

TEST_CASE("Read/Write uint32_t values") {
    std::vector<uint32_t> source_data = {0x44332211, 0x11223344, 0, 1, 2, 3, 4};
    std::vector<uint32_t> destination_data;
    std::vector<char> stream_data(sizeof(uint32_t) * source_data.size());
    MockStream<char> stream(stream_data);
    // Write data to stream.
    for (const auto &value : source_data) {
        CHECK(Serialization::write_uint32(stream, value));
    }
    // Read data from stream.
    for (const auto &value : source_data) {
        uint32_t read_value;
        CHECK(Serialization::read_uint32(stream, &read_value));
        CHECK(read_value == value);
    }
}

TEST_CASE("Read/Write uint16_t values") {
    std::vector<uint16_t> source_data = {0x2211, 0x1122, 0, 1, 2, 3, 4};
    std::vector<uint16_t> destination_data;
    std::vector<char> stream_data(sizeof(uint16_t) * source_data.size());
    MockStream<char> stream(stream_data);
    // Write data to stream.
    for (const auto &value : source_data) {
        CHECK(Serialization::write_uint16(stream, value));
    }
    // Read data from stream.
    for (const auto &value : source_data) {
        uint16_t read_value;
        CHECK(Serialization::read_uint16(stream, &read_value));
        CHECK(read_value == value);
    }
}

TEST_CASE("Read/Write uint8_t values") {
    std::vector<uint8_t> source_data = {0x00, 0xFF, 0, 1, 2, 3, 4};
    std::vector<uint8_t> destination_data;
    std::vector<char> stream_data(sizeof(uint8_t) * source_data.size());
    MockStream<char> stream(stream_data);
    // Write data to stream.
    for (const auto &value : source_data) {
        CHECK(Serialization::write_uint8(stream, value));
    }
    // Read data from stream.
    for (const auto &value : source_data) {
        uint8_t read_value;
        CHECK(Serialization::read_uint8(stream, &read_value));
        CHECK(read_value == value);
    }
}
