#include "centroid/centroid_files.hpp"

bool Centroid::Files::Bpks::write_header(std::ostream &stream,
                                        const Header &header) {
    //auto footer_size = static_cast<char>(sizeof(Grid::Parameters) +
                                         //sizeof(Grid::Files::Dat::Parameters));
    //Grid::Files::Dat::Parameters file_parameters = {1, footer_size};
    //stream.write(reinterpret_cast<const char *>(&parameters),
                 //sizeof(parameters));
    //stream.write(reinterpret_cast<const char *>(&file_parameters),
                 //sizeof(file_parameters));
    return stream.good();
}
