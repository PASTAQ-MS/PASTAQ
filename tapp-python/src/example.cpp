#include <iostream>

#include "grid/grid.hpp"
#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

namespace py = pybind11;

struct Chromatogram {
    std::vector<double> intensity;
    std::vector<double> rt;
};

struct RawData {
    std::vector<Grid::Point> data;
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
        std::stable_sort(points.begin(), points.end(), sort_points);

        // Find the tic.
        Chromatogram tic = {};
        double previous_rt = -std::numeric_limits<double>::infinity();
        int i = -1;
        for (const auto &point : points) {
            if (point.rt > previous_rt || i == -1) {
                tic.intensity.push_back(0);
                tic.rt.push_back(point.rt);
                ++i;
            }
            tic.intensity[i] += point.value;
        }
        return tic;
    }
};

RawData raw_data(
    // TODO: string arguments for the file name?
) {
    return RawData{
        {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {0, 1, 2}},
    };
}

PYBIND11_MODULE(example, m) {
    // Documentation.
    m.doc() = "pybind11 example plugin";

    // Structs.
    py::class_<Grid::Point>(m, "RawPoint")
        .def_readwrite("mz", &Grid::Point::mz)
        .def_readwrite("rt", &Grid::Point::rt)
        .def_readwrite("value", &Grid::Point::value);

    py::class_<Chromatogram>(m, "Chromatogram")
        .def_readwrite("intensity", &Chromatogram::intensity)
        .def_readwrite("rt", &Chromatogram::rt);

    py::class_<RawData>(m, "RawData", py::buffer_protocol())
        .def_readwrite("data", &RawData::data)
        .def("tic", &RawData::tic)
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
        });

    // Functions.
    m.def("raw_data", &raw_data, "Read raw data from the given mzXML file");
}
