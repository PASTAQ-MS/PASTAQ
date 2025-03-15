#include "timsdata_pybind.h"
#include "tims_visualization.h"
#include <pybind11/stl.h>
#include <pybind11/numpy.h> 

using namespace timsdata;
using namespace timsvis;


void register_tims_data(py::module_ &m) {
    py::enum_<pressure_compensation_strategy>(m, "PressureCompensationStrategy")
        .value("NoPressureCompensation", pressure_compensation_strategy::NoPressureCompensation)
        .value("AnalyisGlobalPressureCompensation", pressure_compensation_strategy::AnalyisGlobalPressureCompensation)
        .value("PerFramePressureCompensation", pressure_compensation_strategy::PerFramePressureCompensation)
        .value("PerFramePressureCompensationWithMissingReference", pressure_compensation_strategy::PerFramePressureCompensationWithMissingReference)
        .export_values();

    // Expose Spectrum Class
    py::class_<TimsSpectrum>(m, "TimsSpectrum")
        .def("getMz", &TimsSpectrum::getMz, "Get m/z values.")
        .def("getIntensity", &TimsSpectrum::getIntensity, "Get intensity values.")
        .def("getMobility", &TimsSpectrum::getMobility, "Get mobility.")
        .def("getNumberOfPeaks", &TimsSpectrum::getNumberOfPeaks, "Get number of peaks.");

    // Expose Scan Class
    // py::class_<TimsScan>(m, "TimsScan")
    //     .def("getSpectra", &TimsScan::getSpectra, "Get list of spectra.")
    //     .def("getNumberOfSpectra", &TimsScan::getNumberOfSpectra, "Get number of spectra.");

    // Expose Frame Class
    py::class_<Frame>(m, "Frame")
        .def("getTime", &Frame::getTime, "Get time.")
        .def("getNumScans", &Frame::getNumScans, "Get number of scans.")
        .def("getSpectra",&Frame::getSpectra,"Get spectra as a list");
    
    py::class_<PyFrameProxy>(m, "PyFrameProxy")
        .def("getNbrScans", &PyFrameProxy::getNbrScans)
        .def("getTotalNbrPeaks", &PyFrameProxy::getTotalNbrPeaks)
        .def("getNbrPeaks", &PyFrameProxy::getNbrPeaks)
        .def("get_data", &PyFrameProxy::get_data, "Returns scan data as a list")
        .def("getScanX", &PyFrameProxy::getScanX)
        .def("getScanY", &PyFrameProxy::getScanY);
    
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
        .def("getTdfFile", &timsdata::TimsData::getTdfFile, "Get the path to the .tdf SQLite file.")
        .def("readScans", [](TimsData &self, int64_t frame_id, uint32_t scan_begin, uint32_t scan_end) {
            FrameProxy result = self.readScans(frame_id, scan_begin, scan_end);
            return PyFrameProxy(result.getNbrScans(), result.getData());  // Convert to Python-friendly structure
        }, "Reads TIMS scan data")
        // Exposing conversion functions
        .def("indexToMz", [](TimsData& self, int64_t frame_id, const std::vector<double>& in) {
            std::vector<double> out;
            self.indexToMz(frame_id, in, out);
            return out;
        }, "Convert index to m/z values.")

        .def("mzToIndex", [](TimsData& self, int64_t frame_id, const std::vector<double>& in) {
            std::vector<double> out;
            self.mzToIndex(frame_id, in, out);
            return out;
        }, "Convert m/z to index values.")

        .def("scanNumToOneOverK0", [](TimsData& self, int64_t frame_id, const std::vector<double>& in) {
            std::vector<double> out;
            self.scanNumToOneOverK0(frame_id, in, out);
            return out;
        }, "Convert scan number to 1/K0.")

        .def("oneOverK0ToScanNum", [](TimsData& self, int64_t frame_id, const std::vector<double>& in) {
            std::vector<double> out;
            self.oneOverK0ToScanNum(frame_id, in, out);
            return out;
        }, "Convert 1/K0 to scan number.")

        .def("scanNumToVoltage", [](TimsData& self, int64_t frame_id, const std::vector<double>& in) {
            std::vector<double> out;
            self.scanNumToVoltage(frame_id, in, out);
            return out;
        }, "Convert scan number to voltage.")

        .def("voltageToScanNum", [](TimsData& self, int64_t frame_id, const std::vector<double>& in) {
            std::vector<double> out;
            self.voltageToScanNum(frame_id, in, out);
            return out;
        }, "Convert voltage to scan number.")
        
        .def("populateFrames", 
             &TimsData::populateFrames, 
             py::arg("filter") = "", 
             "Populate frames, optionally filtering by an SQL condition.")

        .def("selectFrames", &TimsData::selectFrames,
            py::arg("min_id") = std::nullopt,
            py::arg("max_id") = std::nullopt,
            py::arg("min_time") = std::nullopt,
            py::arg("max_time") = std::nullopt,
            py::arg("polarity") = std::nullopt,
            py::arg("msms_type") = std::nullopt,
            py::arg("scan_mode") = std::nullopt,
            py::arg("tims_id") = std::nullopt,
            py::arg("min_max_intensity") = std::nullopt,
            py::arg("max_max_intensity") = std::nullopt,
            py::arg("min_summed_intensities") = std::nullopt,
            py::arg("max_summed_intensities") = std::nullopt,
            py::arg("min_num_peaks") = std::nullopt,
            py::arg("max_num_peaks") = std::nullopt)

        .def("populateSpectra", 
            [](timsdata::TimsData& self, 
                py::object min_mz, py::object max_mz, 
                py::object min_mobility, py::object max_mobility) {
         
                     // Convert Python None to std::optional<double>
                     std::optional<double> min_mz_opt = min_mz.is_none() ? std::nullopt : std::make_optional(min_mz.cast<double>());
                     std::optional<double> max_mz_opt = max_mz.is_none() ? std::nullopt : std::make_optional(max_mz.cast<double>());
                     std::optional<double> min_mobility_opt = min_mobility.is_none() ? std::nullopt : std::make_optional(min_mobility.cast<double>());
                     std::optional<double> max_mobility_opt = max_mobility.is_none() ? std::nullopt : std::make_optional(max_mobility.cast<double>());

                 self.populateSpectra(min_mz_opt, max_mz_opt, min_mobility_opt, max_mobility_opt);
                 },
             py::arg("min_mz") = py::none(),
             py::arg("max_mz") = py::none(),
             py::arg("min_mobility") = py::none(),
             py::arg("max_mobility") = py::none(),
             "Populate spectra for the previously loaded frames, with optional filtering by m/z and mobility.")
        
        .def("getFrames", &timsdata::TimsData::getFrames, "Get all frames.");

        py::class_<timsvis::HeatmapData>(m, "HeatmapData")
        .def_readonly("x_values", &timsvis::HeatmapData::x_edges, 
                  "List of x-values (m/z) corresponding to each column in the heatmap.")
        .def_readonly("y_values", &timsvis::HeatmapData::y_edges, 
                  "List of y-values (Time) corresponding to each row in the heatmap.")
        .def_readonly("z_values", &timsvis::HeatmapData::value_matrix, 
                  "List of intensity values (log-transformed) at each (x, y) position in the heatmap.");

        py::class_<timsvis::TimsVisualization>(m, "TimsVisualization")
        .def(py::init<const timsdata::TimsData&>(), 
         "Constructor for TimsVisualization class, requires a TimsData object to initialize.")
        .def("heatmap_data", [](timsvis::TimsVisualization& self, int x_bins, int y_bins) {
            auto result = self.heatmap_data(x_bins, y_bins);
            return py::make_tuple(
            py::array_t<double>(result.x_edges.size(), result.x_edges.data()),  // x_values
            py::array_t<double>(result.y_edges.size(), result.y_edges.data()),  // y_values
            py::array_t<double>({result.value_matrix.size(), result.value_matrix[0].size()}, &result.value_matrix[0][0])   // z_values
            );
        }, 
        "Generates heatmap data with x_bins and y_bins, returning x_values, y_values, and z_values as NumPy arrays.");
    
    // py::class_<timsvis::HeatmapData>(m, "HeatmapData")
    //     .def_readonly("value_matrix", &HeatmapData::value_matrix)
    //     .def_readonly("x_edges", &HeatmapData::x_edges)
    //     .def_readonly("y_edges", &HeatmapData::y_edges);

    // py::class_<timsvis::TimsVisualization>(m, "TimsVisualization")
    //     .def(py::init<const timsdata::TimsData&>())
    //     .def("heatmap_data", [](timsvis::TimsVisualization& self, int x_bins, int y_bins) {
    //         auto result = self.heatmap_data(x_bins, y_bins);
    //         return py::make_tuple(
    //             py::array_t<double>({result.value_matrix.size(), result.value_matrix[0].size()}, &result.value_matrix[0][0]),
    //             py::array_t<double>(result.x_edges.size(), result.x_edges.data()),
    //             py::array_t<double>(result.y_edges.size(), result.y_edges.data())
    //         );
    //     });
}
