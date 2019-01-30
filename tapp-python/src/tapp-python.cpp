#include <filesystem>
#include <fstream>
#include <iostream>

#include "grid/grid.hpp"
#include "grid/grid_files.hpp"
#include "grid/xml_reader.hpp"
#include "centroid/centroid.hpp"
#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

namespace py = pybind11;

struct Chromatogram {
    std::vector<double> retention_time;
    std::vector<double> intensity;
};


std::vector<std::tuple<int,int,int>> find_local_maxima(uint64_t n, uint64_t m, std::vector<double> &data) {
    std::vector<std::tuple<int,int,int>> ret;
    auto local_max = Centroid::find_local_maxima({0,0, {{n,m},{},{}, Instrument::QUAD, 0x00}}, data);
    for (const auto &p : local_max) {
        std::tuple<int,int,int> point;
        std::get<0>(point) = p.i;
        std::get<1>(point) = p.j;
        std::get<2>(point) = p.value;
        ret.push_back(point);
    }
    return ret;
}

std::vector<std::tuple<int,int,double,double,double,double,double,double,double>> find_peaks(
        uint64_t n,
        uint64_t m,
        double min_mz,
        double max_mz,
        double min_rt,
        double max_rt,
        std::vector<double> &data
        ) {
    std::vector<std::tuple<int,int,double,double,double,double,double,double,double>> ret;
    Centroid::Parameters parameters = {0,0, {{n,m},{min_rt, max_rt, min_mz, max_mz},{}, Instrument::ORBITRAP, 0x00}};
    auto local_max = Centroid::find_local_maxima(parameters, data);
    for (const auto &p : local_max) {
        auto peak = Centroid::build_peak(p, parameters, data);
        if (peak) {
            std::tuple<int,int,double,double,double,double,double,double,double> x;
            std::get<0>(x) = peak.value().i;
            std::get<1>(x) = peak.value().j;
            std::get<2>(x) = peak.value().mz;
            std::get<3>(x) = peak.value().rt;
            std::get<4>(x) = peak.value().height;
            std::get<5>(x) = peak.value().total_intensity;
            std::get<6>(x) = peak.value().sigma_mz;
            std::get<7>(x) = peak.value().sigma_rt;
            std::get<8>(x) = peak.value().border_background;
            ret.push_back(x);
        }
    }
    return ret;
}

struct RawData {
    std::vector<Grid::Point> data;
    std::string file_name;
    // TODO: Store dimensions, instrument, min/max mz/rt

    // Extract the total ion chromatogram from the RawData.
    Chromatogram tic() {
        // Make a copy of the original data before sorting.
        std::vector<Grid::Point> points = data;

        // Sort the points by retention time.
        auto sort_points = [](const Grid::Point &p1,
                              const Grid::Point &p2) -> bool {
            return (p1.rt < p2.rt);
        };
        std::sort(points.begin(), points.end(), sort_points);

        // Find the tic.
        Chromatogram tic = {};
        double previous_rt = -std::numeric_limits<double>::infinity();
        int i = -1;
        for (const auto &point : points) {
            if (point.rt > previous_rt || i == -1) {
                previous_rt = point.rt;
                tic.intensity.push_back(0);
                tic.retention_time.push_back(point.rt);
                ++i;
            }
            tic.intensity[i] += point.value;
        }
        return tic;
    }

    // TODO(alex): Save RawData object to binary format rawdump.
    void dump(std::string file_name) {
        // Open file stream.
        std::filesystem::path output_file = file_name;
        std::ofstream stream;
        stream.open(output_file);
        if (!stream) {
            std::ostringstream error_stream;
            error_stream << "error: couldn't open output file" << output_file;
            throw std::invalid_argument(error_stream.str());
        }

        if (!Grid::Files::Rawdump::write(stream, data)) {
            std::ostringstream error_stream;
            error_stream << "error: couldn't write file succesfully"
                         << output_file;
            throw std::invalid_argument(error_stream.str());
        }
    }
    void load(std::string file_name) {
        // Open file stream.
        std::filesystem::path input_file = file_name;
        std::ifstream stream;
        stream.open(input_file);
        if (!stream) {
            std::ostringstream error_stream;
            error_stream << "error: couldn't open input file" << input_file;
            throw std::invalid_argument(error_stream.str());
        }

        // TODO(alex): load object from binary format rawdump
        if (!Grid::Files::Rawdump::read(stream, data)) {
            data = {};
        }

        this->file_name = file_name;
    }
    std::vector<double> mz() {
        std::vector<double> ret = std::vector<double>(data.size());
        for (size_t i = 0; i < data.size(); ++i) {
            ret[i] = data[i].mz;
        }
        return ret;
    }
    std::vector<double> rt() {
        std::vector<double> ret = std::vector<double>(data.size());
        for (size_t i = 0; i < data.size(); ++i) {
            ret[i] = data[i].rt;
        }
        return ret;
    }
    std::vector<double> intensity() {
        std::vector<double> ret = std::vector<double>(data.size());
        for (size_t i = 0; i < data.size(); ++i) {
            ret[i] = data[i].value;
        }
        return ret;
    }
};

//std::vector<double> splat_points(
    //std::vector<Grid::Point> points,
    //double min_mz,
    //double max_mz,
    //double min_rt,
    //double max_mz,
    //double delta_rt,
    //double delta_mz,
    //std::tuple<double, double> sigma_mz,
    //double sigma_rt
    //) {
//}

// IMPORTANT: This is awful design...
RawData raw_data_load_dump(std::string file_name) {
    // TODO(alex): Check for proper file extension.
    RawData raw_data = {};
    raw_data.load(file_name);
    return raw_data;
}

// TODO(alex): Should this be encapsulated in a RawData or raw_data namespace?
RawData raw_data(std::string file_name, double min_mz, double max_mz,
                 double min_rt, double max_rt) {
    // Setup infinite range if no point was specified.
    min_rt = min_rt < 0 ? 0 : min_rt;
    max_rt = max_rt < 0 ? std::numeric_limits<double>::infinity() : max_rt;
    min_mz = min_mz < 0 ? 0 : min_mz;
    max_mz = max_mz < 0 ? std::numeric_limits<double>::infinity() : max_mz;

    // Sanity check the min/max rt/mz.
    if (min_rt >= max_rt) {
        std::ostringstream error_stream;
        error_stream << "error: min_rt >= max_rt (min_rt: " << min_rt
                     << ", max_rt: " << max_rt << ")";
        throw std::invalid_argument(error_stream.str());
    }
    if (min_mz >= max_mz) {
        std::ostringstream error_stream;
        error_stream << "error: min_mz >= max_mz (min_mz: " << min_mz
                     << ", max_mz: " << max_mz << ")";
        throw std::invalid_argument(error_stream.str());
    }

    // TODO(alex): Check for proper file extension.

    std::filesystem::path input_file = file_name;

    // Open file stream.
    std::ifstream stream;
    stream.open(input_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }

    // Read the data from the file.
    Grid::Parameters parameters = {
        {}, {min_rt, max_rt, min_mz, max_mz}, {}, 0x00, 0x00,
    };
    auto points = XmlReader::read_next_scan(stream, parameters);
    std::vector<Grid::Point> all_points = {};
    do {
        all_points.insert(end(all_points), begin(points.value()),
                          end(points.value()));
        points = XmlReader::read_next_scan(stream, parameters);
    } while (points != std::nullopt);

    return RawData{
        all_points,
        input_file.filename(),
    };
}

PYBIND11_MODULE(tapp, m) {
    // Documentation.
    m.doc() = "tapp documentation";

    // Structs.
    py::class_<Grid::Point>(m, "RawPoint")
        .def_readonly("mz", &Grid::Point::mz)
        .def_readonly("rt", &Grid::Point::rt)
        .def_readonly("value", &Grid::Point::value)
        .def("__repr__", [](const Grid::Point &p) {
            return "RawPoint {mz: " + std::to_string(p.mz) +
                   ", rt: " + std::to_string(p.rt) +
                   ", intensity: " + std::to_string(p.value) + "}";
        });

    py::class_<Chromatogram>(m, "Chromatogram")
        .def_readonly("intensity", &Chromatogram::intensity)
        .def_readonly("retention_time", &Chromatogram::retention_time);

    py::class_<RawData>(m, "RawData", py::buffer_protocol())
        .def_readonly("data", &RawData::data)
        .def("tic", &RawData::tic)
        .def("dump", &RawData::dump)
        .def("load", &RawData::load)
        .def("mz", &RawData::mz)
        .def("rt", &RawData::rt)
        .def("intensity", &RawData::intensity)
        .def_buffer([](RawData &m) -> py::buffer_info {
            return py::buffer_info(
                // Pointer to buffer.
                m.data.data(),
                // Size of one scalar.
                sizeof(double),
                // Format descriptor.
                py::format_descriptor<double>::format(),
                // Number of dimensions.
                2,
                // Number of elements in each dimension.
                {int(m.data.size()), 3},
                // Stride for each dimension (bytes).
                {sizeof(double) * (m.data.size() - 1), sizeof(double)});
        })
        .def("__repr__", [](const RawData &rd) {
            return "RawData {file_name: " + rd.file_name +
                   ", number_of_points: " + std::to_string(rd.data.size()) +
                   "}";
        });

    // Functions.
    m.def("find_local_maxima", &find_local_maxima, "Find local maxima on the given smoothed grid",
          py::arg("num_mz"),
          py::arg("num_rt"),
          py::arg("grid"));
    m.def("find_peaks", &find_peaks, "Find peaks on the given smoothed grid",
          py::arg("num_mz"),
          py::arg("num_rt"),
          py::arg("min_mz"),
          py::arg("max_mz"),
          py::arg("min_rt"),
          py::arg("max_rt"),
          py::arg("grid"));
    m.def("raw_data", &raw_data, "Read raw data from the given mzXML file",
          py::arg("file_name"), py::arg("min_mz") = -1.0,
          py::arg("max_mz") = -1.0, py::arg("min_rt") = -1.0,
          py::arg("max_rt") = -1.0);
    m.def("raw_data_load_dump", &raw_data_load_dump,
          "Read raw data from the given rawdump file", py::arg("file_name"));
    m.def("dummy_test",
          []() {
              return RawData{{{1, 1.0, 1}, {1, 1.0, 1}, {1, 2.0, 2}},
                             "dummy.mzXML"};
          },
          "DEBUG: dummy test");
}
