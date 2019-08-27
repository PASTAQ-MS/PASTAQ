#ifndef GRID_RAWDATASERIALIZE_HPP
#define GRID_RAWDATASERIALIZE_HPP

#include <iostream>

#include "raw_data.hpp"

// This namespace groups the functions used to serialize RawData data structures
// into a binary stream.
namespace RawData::Serialize {

// RawData::Scan
bool read_scan(std::istream &stream, Scan *scan);
bool write_scan(std::ostream &stream, const Scan &scan);

// RawData::PrecursorInformation
bool read_precursor_info(std::istream &stream,
                         PrecursorInformation *precursor_info);
bool write_precursor_info(std::ostream &stream,
                          const PrecursorInformation &precursor_info);

// RawData::RawData
bool read_raw_data(std::istream &stream, RawData *raw_data);
bool write_raw_data(std::ostream &stream, const RawData &raw_data);

}  // namespace RawData::Serialize

#endif /* GRID_RAWDATASERIALIZE_HPP */
