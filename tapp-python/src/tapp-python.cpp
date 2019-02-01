#include <filesystem>
#include <fstream>
#include <iostream>

#include "centroid/centroid.hpp"
#include "grid/grid.hpp"
#include "grid/grid_files.hpp"
#include "grid/xml_reader.hpp"
#include "pybind11/numpy.h"
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

namespace py = pybind11;

namespace PythonAPI {
RawData::RawData read_mzxml(std::string file_name, double min_mz, double max_mz,
                            double min_rt, double max_rt,
                            std::string instrument_type_str,
                            double resolution_ms1, double resolution_msn,
                            double reference_mz, std::string polarity_str) {
    // Setup infinite range if no point was specified.
    min_rt = min_rt < 0 ? 0 : min_rt;
    max_rt = max_rt < 0 ? std::numeric_limits<double>::infinity() : max_rt;
    min_mz = min_mz < 0 ? 0 : min_mz;
    max_mz = max_mz < 0 ? std::numeric_limits<double>::infinity() : max_mz;

    // Parse the instrument type.
    auto instrument_type = Instrument::Type::UNKNOWN;
    // for (auto &ch : instrument_type_str) {
    // ch = std::tolower(ch);
    //}
    if (instrument_type_str == "orbitrap") {
        instrument_type = Instrument::Type::ORBITRAP;
    } else if (instrument_type_str == "tof") {
        instrument_type = Instrument::Type::TOF;
    } else if (instrument_type_str == "quad") {
        instrument_type = Instrument::Type::QUAD;
    } else if (instrument_type_str == "fticr") {
        instrument_type = Instrument::Type::FTICR;
    } else {
        std::ostringstream error_stream;
        error_stream << "the given instrument is not supported";
        throw std::invalid_argument(error_stream.str());
    }
    // Parse the polarity.
    auto polarity = RawData::Polarity::BOTH;
    // for (auto &ch : polarity_str) {
    // ch = std::tolower(ch);
    //}
    if (polarity_str == "" || polarity_str == "both" || polarity_str == "+-" ||
        polarity_str == "-+") {
        polarity = RawData::Polarity::BOTH;
    } else if (polarity_str == "+" || polarity_str == "pos" ||
               polarity_str == "positive") {
        polarity = RawData::Polarity::POSITIVE;
    } else if (polarity_str == "-" || polarity_str == "neg" ||
               polarity_str == "negative") {
        polarity = RawData::Polarity::NEGATIVE;
    } else {
        std::ostringstream error_stream;
        error_stream << "the given polarity is not supported. choose "
                        "between '+', '-', 'both' (default)";
        throw std::invalid_argument(error_stream.str());
    }

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

    std::filesystem::path input_file = file_name;

    // Check for proper file extension.
    std::string extension = input_file.extension();
    std::string lowercase_extension = extension;
    for (auto &ch : lowercase_extension) {
        ch = std::tolower(ch);
    }
    if (lowercase_extension != ".mzxml") {
        std::ostringstream error_stream;
        error_stream << "invalid file type: expected 'mzXML' but given '"
                     << extension << "'";
        throw std::invalid_argument(error_stream.str());
    }

    // Open file stream.
    std::ifstream stream;
    stream.open(input_file);
    if (!stream) {
        std::ostringstream error_stream;
        error_stream << "error: couldn't open input file" << input_file;
        throw std::invalid_argument(error_stream.str());
    }

    auto raw_data = XmlReader::read_mzxml(
        stream, min_mz, max_mz, min_rt, max_rt, instrument_type, resolution_ms1,
        resolution_msn, reference_mz, polarity);
    if (!raw_data) {
        std::ostringstream error_stream;
        error_stream << "error: an error occurred when reading the file"
                     << input_file;
        throw std::invalid_argument(error_stream.str());
    }
    return raw_data.value();
}

std::string to_string(const Instrument::Type &instrument_type) {
    switch (instrument_type) {
        case Instrument::Type::QUAD:
            return "QUAD";
        case Instrument::Type::TOF:
            return "TOF";
        case Instrument::Type::FTICR:
            return "FTICR";
        case Instrument::Type::ORBITRAP:
            return "ORBITRAP";
        case Instrument::Type::UNKNOWN:
            return "UNKNOWN";
    };
    return "UNKNOWN";
}

std::string to_string(const RawData::Polarity &polarity) {
    switch (polarity) {
        case RawData::Polarity::POSITIVE:
            return "POSITIVE";
        case RawData::Polarity::NEGATIVE:
            return "NEGATIVE";
        case RawData::Polarity::BOTH:
            return "BOTH";
    };
    return "UNKNOWN";
}

}  // namespace PythonAPI

PYBIND11_MODULE(tapp, m) {
    // Documentation.
    m.doc() = "tapp documentation";

    // Structs.
    py::class_<RawData::Scan>(m, "Scan")
        .def_readonly("scan_number", &RawData::Scan::scan_number)
        .def_readonly("ms_level", &RawData::Scan::ms_level)
        .def_readonly("num_points", &RawData::Scan::num_points)
        .def_readonly("retention_time", &RawData::Scan::retention_time)
        .def_readonly("mz", &RawData::Scan::mz)
        .def_readonly("intensity", &RawData::Scan::intensity)
        .def_readonly("polarity", &RawData::Scan::polarity);

    py::class_<Instrument::Type>(m, "Instrument")
        .def("__repr__", [](const Instrument::Type &instrument_type) {
            return PythonAPI::to_string(instrument_type);
        });

    py::class_<RawData::Polarity>(m, "Polarity")
        .def("__repr__", [](const RawData::Polarity &polarity) {
            return PythonAPI::to_string(polarity);
        });

    py::class_<RawData::RawData>(m, "RawData")
        .def_readonly("scans", &RawData::RawData::scans)
        .def_readonly("instrument_type", &RawData::RawData::instrument_type)
        .def_readonly("resolution_ms1", &RawData::RawData::resolution_ms1)
        .def_readonly("resolution_msn", &RawData::RawData::resolution_msn)
        .def_readonly("reference_mz", &RawData::RawData::reference_mz)
        .def_readonly("min_mz", &RawData::RawData::min_mz)
        .def_readonly("max_mz", &RawData::RawData::max_mz)
        .def_readonly("min_rt", &RawData::RawData::min_rt)
        .def_readonly("max_rt", &RawData::RawData::max_rt)
        .def("__repr__", [](const RawData::RawData &rd) {
            return "RawData:\n> instrument_type: " +
                   PythonAPI::to_string(rd.instrument_type) +
                   "\n> resolution_ms1: " + std::to_string(rd.resolution_ms1) +
                   "\n> resolution_msn: " + std::to_string(rd.resolution_msn) +
                   "\n> reference_mz: " + std::to_string(rd.reference_mz) +
                   "\n> min_mz: " + std::to_string(rd.min_mz) +
                   "\n> max_mz: " + std::to_string(rd.max_mz) +
                   "\n> min_rt: " + std::to_string(rd.min_rt) +
                   "\n> max_rt: " + std::to_string(rd.max_rt) +
                   "\n> number of scans: " + std::to_string(rd.scans.size());
        });

    // Functions.
    m.def("read_mzxml", &PythonAPI::read_mzxml,
          "Read raw data from the given mzXML file ", py::arg("file_name"),
          py::arg("min_mz") = -1.0, py::arg("max_mz") = -1.0,
          py::arg("min_rt") = -1.0, py::arg("max_rt") = -1.0,
          py::arg("instrument_type") = "", py::arg("resolution_ms1"),
          py::arg("resolution_msn"), py::arg("reference_mz"),
          py::arg("polarity") = "");
}
