#include <iostream>
#include <vector>

#include "tapp-grid.hpp"

// TODO: This binary file will contain the data at the beginning followed by
// a header with the parameters that we need to interpret the data properly.
// Hopefully doing it this way will make it so it is compatible with the
// existing solutions.
// TODO: This reader/writer module should deal properly with endianness. We will
// save the data as little-endian and if we encounter an architecture that is
// big endian we will swap the bytes around.
namespace DatFile {

// Load an entire file from the binary stream into the destination vector and
// the parameters structure.
bool load(std::istream &stream, std::vector<double> &destination,
          Grid::Parameters &parameters);

// Load the parameters from the footer of the binary straem.
std::optional<Grid::Parameters> load_parameters(std::istream &stream);

// Write the entire source vector and parameters struct into the destination
// binary stream.
bool write(std::ostream &stream, const std::vector<double> &source,
           const Grid::Parameters &parameters);

// TODO(alex): Allow for reading or loading specific ranges.
// bool load_range(std::istream &stream, const Grid::Bounds &bounds,
// std::vector<double> &destination, Grid::Parameters &parameters);
// bool write_range(std::ostream &stream, const Grid::Bounds &bounds,
// const std::vector<double> &source,
// const Grid::Parameters &parameters);

// Read/Write an unsigned integer from/to the stream.
bool load_uint32(std::istream &stream, uint32_t *i);
bool save_uint32(std::ostream &stream, uint32_t i);

}  // namespace DatFile
