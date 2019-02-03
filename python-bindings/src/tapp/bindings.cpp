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
    for (auto &ch : instrument_type_str) {
        ch = std::tolower(ch);
    }
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
    for (auto &ch : polarity_str) {
        ch = std::tolower(ch);
    }
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

std::tuple<uint64_t, uint64_t> calculate_dimensions(
    const RawData::RawData &raw_data, double avg_rt_fwhm,
    uint64_t num_samples_per_peak_mz, uint64_t num_samples_per_peak_rt) {
    // Calculate the number of sampling points in the rt dimension.
    //
    // NOTE(alex): Since the average retention time is given in FWHM and under
    // the assumption of Gaussian chromatographic peaks, the FWHM ≈ 2.355 *
    // sigma. We need then 3 sigma left and right of the center of the
    // chromatographic peak to cover a 99.7 % of the gaussian peak area.
    double sigma_rt = avg_rt_fwhm / (2 * std::sqrt(2 * std::log(2)));
    double base_width_rt = sigma_rt * 6;
    double delta_rt = base_width_rt / num_samples_per_peak_rt;
    uint64_t num_points_rt =
        std::ceil((raw_data.max_rt - raw_data.min_rt) / delta_rt);

    double fwhm_ref = raw_data.reference_mz / raw_data.resolution_ms1;

    // FIXME: This only works for ORBITRAP data for now.
    uint64_t num_points_mz =
        num_samples_per_peak_mz * 2 * std::pow(raw_data.reference_mz, 1.5) /
        fwhm_ref *
        (1 / std::sqrt(raw_data.min_mz) - 1 / std::sqrt(raw_data.max_mz));
    return std::tuple<uint64_t, uint64_t>(num_points_mz + 1, num_points_rt + 1);
}

double mz_at(const RawData::RawData &raw_data, uint64_t num_samples_per_peak_mz,
             uint64_t n) {
    // FIXME: This only works for ORBITRAP data for now.
    double a = 1 / std::sqrt(raw_data.min_mz);
    double fwhm_ref = raw_data.reference_mz / raw_data.resolution_ms1;
    double b = fwhm_ref / std::pow(raw_data.reference_mz, 1.5) * n / 2 /
               num_samples_per_peak_mz;
    double c = a - b;
    return 1 / (c * c);
}

double fwhm_at(const RawData::RawData &raw_data, double mz) {
    // FIXME: This only works for ORBITRAP data for now.
    double fwhm_ref = raw_data.reference_mz / raw_data.resolution_ms1;
    return fwhm_ref * std::pow(mz / raw_data.reference_mz, 1.5);
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

struct Mesh {
    uint64_t n;
    uint64_t m;
    std::vector<double> matrix;
    std::vector<double> bins_mz;
    std::vector<double> bins_rt;
};

Mesh resample(const RawData::RawData &raw_data, double avg_rt_fwhm,
              uint64_t num_samples_per_peak_mz,
              uint64_t num_samples_per_peak_rt) {
    auto [n, m] =
        calculate_dimensions(raw_data, avg_rt_fwhm, num_samples_per_peak_mz,
                             num_samples_per_peak_rt);
    Mesh mesh;
    mesh.n = n;
    mesh.m = m;
    mesh.matrix = std::vector<double>(n * m);
    mesh.bins_mz = std::vector<double>(n);
    mesh.bins_rt = std::vector<double>(m);

    // Generate bins_mz.
    for (size_t i = 0; i < n; ++i) {
        mesh.bins_mz[i] = mz_at(raw_data, num_samples_per_peak_mz, i);
    }
    // Generate bins_rt.
    double delta_rt = (raw_data.max_rt - raw_data.min_rt) / (m - 1);
    for (size_t j = 0; j < m; ++j) {
        mesh.bins_rt[j] = raw_data.min_rt + delta_rt * j;
    }

    double sigma_rt = avg_rt_fwhm / 2.355; // FIXME: Approx

    // DEBUG
    std::cout << "raw_data.min_mz: " << raw_data.min_mz << std::endl;
    std::cout << "raw_data.max_mz: " << raw_data.max_mz << std::endl;
    std::cout << "mesh.bins_mz[0]: " << mesh.bins_mz[0] << std::endl;
    std::cout << "mesh.bins_mz[n-1]: " << mesh.bins_mz[n - 1] << std::endl;
    std::cout << "raw_data.min_rt: " << raw_data.min_rt << std::endl;
    std::cout << "raw_data.max_rt: " << raw_data.max_rt << std::endl;
    std::cout << "mesh.bins_rt[0]: " << mesh.bins_rt[0] << std::endl;
    std::cout << "mesh.bins_rt[m-1]: " << mesh.bins_rt[m - 1] << std::endl;

    for (const auto &scan : raw_data.scans) {
        // double sigma_rt = avg_rt_fwhm / 2.355;
        // Calculate the min and max indexes for retention time.

        // Find the bin for the current retention time.
        // NOTE: y_index.
        double current_rt = scan.retention_time;
        size_t index_rt = (current_rt - raw_data.min_rt) / delta_rt;

        // The smoothing kernel in rt is +-(num_samples_per_peak_rt/2).
        int64_t j_min = index_rt - num_samples_per_peak_rt / 2;
        if (j_min < 0) {
            j_min = 0;
        }
        int64_t j_max = index_rt + num_samples_per_peak_rt / 2;
        if (j_max >= m) {
            j_max = m - 1;
        }

        //// DEBUG
        // std::cout << "current_rt: " << rt << std::endl;
        // std::cout << "index_rt: " << index_rt << std::endl;
        // std::cout << "mesh.bins_rt[index_rt]: " << mesh.bins_rt[index_rt]
        //<< std::endl;
        // std::cout << "j_min: " << j_min << std::endl;
        // std::cout << "j_max: " << j_max << std::endl;
        //

        for (size_t k = 0; k < scan.num_points; ++k) {
            double current_intensity = scan.intensity[k];

            // Find the bin for the current mz.
            double current_mz = scan.mz[k];

            // NOTE: x_index
            // FIXME: This only works for ORBITRAP data for now.
            double fwhm_ref = raw_data.reference_mz / raw_data.resolution_ms1;
            uint64_t index_mz =
                num_samples_per_peak_mz * 2 *
                std::pow(raw_data.reference_mz, 1.5) / fwhm_ref *
                (1 / std::sqrt(raw_data.min_mz) - 1 / std::sqrt(current_mz));
            int64_t i_min = index_mz - num_samples_per_peak_mz / 2;
            if (i_min < 0) {
                i_min = 0;
            }
            int64_t i_max = index_mz + num_samples_per_peak_mz / 2;
            if (i_max >= n) {
                i_max = n - 1;
            }
            //// DEBUG
            // std::cout << "current_mz: " << current_mz << std::endl;
            // std::cout << "index_mz: " << index_mz << std::endl;
            // std::cout << "mesh.bins_mz[index_mz]: " << mesh.bins_mz[index_mz]
            //<< std::endl;
            // std::cout << "i_min: " << i_min << std::endl;
            // std::cout << "i_max: " << i_max << std::endl;
            for (size_t j = j_min; j <= j_max; ++j) {
                for (size_t i = i_min; i <= i_max; ++i) {
                    // FIXME: ORBITRAP
                    // NOTE: Should we precalculate this?
                    double sigma_mz = (fwhm_ref * std::pow(current_mz/raw_data.reference_mz, 1.5)) / 2.355; // FIXME: Approx

                    // No need to do boundary check, since we are sure we are
                    // inside the grid.
                    double x = mesh.bins_mz[i];
                    double y = mesh.bins_rt[j];

                    // Calculate the gaussian weight for this point.
                    double a = (x - current_mz) / sigma_mz;
                    double b = (y - current_rt) / sigma_rt;
                    double weight = std::exp(-0.5 * (a * a + b * b));

                    // Set the value, weight and counts.
                    mesh.matrix[i + j * n] += weight * current_intensity;
                }
                // break;
            }
        }
        //break;
    }
    return mesh;
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

    py::class_<PythonAPI::Mesh>(m, "Mesh")
        .def_readonly("n", &PythonAPI::Mesh::n)
        .def_readonly("m", &PythonAPI::Mesh::m)
        .def_readonly("matrix", &PythonAPI::Mesh::matrix)
        .def_readonly("bins_mz", &PythonAPI::Mesh::bins_mz)
        .def_readonly("bins_rt", &PythonAPI::Mesh::bins_rt);

    // Functions.
    m.def("read_mzxml", &PythonAPI::read_mzxml,
          "Read raw data from the given mzXML file ", py::arg("file_name"),
          py::arg("min_mz") = -1.0, py::arg("max_mz") = -1.0,
          py::arg("min_rt") = -1.0, py::arg("max_rt") = -1.0,
          py::arg("instrument_type") = "", py::arg("resolution_ms1"),
          py::arg("resolution_msn"), py::arg("reference_mz"),
          py::arg("polarity") = "")
        .def("calculate_dimensions", &PythonAPI::calculate_dimensions,
             "Calculate the grid parameters for the given raw file",
             py::arg("raw_data"), py::arg("rt_fwhm"), py::arg("num_mz") = 10,
             py::arg("num_rt") = 10)
        .def("mz_at", &PythonAPI::mz_at,
             "Calculate the mz at the given N for the given raw file",
             py::arg("raw_data"), py::arg("num_mz") = 10, py::arg("n"))
        .def("fwhm_at", &PythonAPI::fwhm_at,
             "Calculate the width of the peak at the given m/z for the given "
             "raw file",
             py::arg("raw_data"), py::arg("mz"))
        .def("resample", &PythonAPI::resample,
             "Resample the raw data into a warped grid", py::arg("raw_data"),
             py::arg("rt_fwhm"), py::arg("num_mz") = 10,
             py::arg("num_rt") = 10);
}
