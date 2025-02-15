#include "timsdata_pybind.h"
#include <pybind11/stl.h>

using namespace timsdata;

void register_tims_data(py::module_ &m) {
    py::enum_<pressure_compensation_strategy>(m, "PressureCompensationStrategy")
        .value("NoPressureCompensation", pressure_compensation_strategy::NoPressureCompensation)
        .value("AnalyisGlobalPressureCompensation", pressure_compensation_strategy::AnalyisGlobalPressureCompensation)
        .value("PerFramePressureCompensation", pressure_compensation_strategy::PerFramePressureCompensation)
        .value("PerFramePressureCompensationWithMissingReference", pressure_compensation_strategy::PerFramePressureCompensationWithMissingReference)
        .export_values();
    
    py::class_<timsdata::TimsData>(m, "TimsData")
        .def(py::init<const std::string&, bool, pressure_compensation_strategy>(),
             py::arg("analysis_directory_name"),
             py::arg("use_recalibration") = false,
             py::arg("pressure_compensation") = pressure_compensation_strategy::AnalyisGlobalPressureCompensation)
        .def("getNumberOfFrames", &timsdata::TimsData::getNumberOfFrames, "Get the number of frames in the analysis.")
        .def("getFramesTable", [](const TimsData& self) {
            auto [column_names, frames_table] = self.getFramesTable();

            py::list py_column_names;
            for (const auto& col : column_names) {
                py_column_names.append(py::str(col));
            }
            
            std::vector<std::vector<py::object>> py_frames_table;
            for (const auto& row : frames_table) {
                std::vector<py::object> py_row;
                for (const auto& col : row) {
                    if (std::holds_alternative<int64_t>(col)) {
                        py_row.push_back(py::int_(std::get<int64_t>(col)));
                    } else if (std::holds_alternative<double>(col)) {
                        py_row.push_back(py::float_(std::get<double>(col)));
                    } else if (std::holds_alternative<std::string>(col)) {
                        py_row.push_back(py::str(std::get<std::string>(col)));
                    }
                }
                py_frames_table.push_back(py_row);
            }

            return py::make_tuple(py_column_names, py_frames_table);
        }, "Get the Frames table as a list of rows, where each row is a list of values.")
        .def("getTdfFile", &timsdata::TimsData::getTdfFile, "Get the path to the .tdf SQLite file.");
}
