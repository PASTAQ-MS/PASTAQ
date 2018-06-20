#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <vector>

#include "grid.hpp"
#include "grid_file.hpp"
#include "xml_reader.hpp"

// Type aliases.
using options_map = std::map<std::string, std::string>;

void print_usage() {
    std::cout << "USAGE: grid [-help] [options] <files>" << std::endl;
}

// Helper functions to check if the given string contains a number.
bool is_unsigned_int(std::string& s) {
    std::regex double_regex("^([[:digit:]]+)$");
    return std::regex_search(s, double_regex);
}
bool is_number(std::string& s) {
    std::regex double_regex("^([[:digit:]]+[\\.]?[[:digit:]]*)$");
    return std::regex_search(s, double_regex);
}

// Helper function to trim the whitespace surrounding a string.
void trim_space(std::string& s) {
    auto not_space = [](int ch) { return !std::isspace(ch); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
    s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
}

bool parse_hdr(const std::filesystem::path& path, options_map& options) {
    std::ifstream stream(path);
    std::string parameter;
    const std::string delimiter = "<==>";
    while (stream.good()) {
        std::getline(stream, parameter);
        auto pos = parameter.find(delimiter);
        auto name = parameter.substr(0, pos);
        parameter.erase(0, pos + delimiter.size());
        trim_space(parameter);

        // Bounds.
        if (name == "ConversionStartMass") {
            if (options.find("-min_mz") == options.end()) {
                options["-min_mz"] = parameter;
            }
        } else if (name == "ConversionEndMass") {
            if (options.find("-max_mz") == options.end()) {
                options["-max_mz"] = parameter;
            }
        } else if (name == "ConversionStartTime") {
            if (options.find("-min_rt") == options.end()) {
                options["-min_rt"] = parameter;
            }
        } else if (name == "ConversionEndTime") {
            if (options.find("-max_rt") == options.end()) {
                options["-max_rt"] = parameter;
            }
        } else if (name == "ConversionMassAtSigma") {
            if (options.find("-smooth_mz") == options.end()) {
                options["-smooth_mz"] = parameter;
            }
        } else if (name == "ConversionSigmaMass") {
            if (options.find("-sigma_mz") == options.end()) {
                options["-sigma_mz"] = parameter;
            }
        } else if (name == "ConversionSigmaTime") {
            if (options.find("-sigma_rt") == options.end()) {
                options["-sigma_rt"] = parameter;
            }
        } else if (name == "ConversionMeanDeltaMass") {
            if (options.find("-delta_mz") == options.end()) {
                options["-delta_mz"] = parameter;
            }
        } else if (name == "ConversionMeanDeltaTime") {
            if (options.find("-delta_rt") == options.end()) {
                options["-delta_rt"] = parameter;
            }
        } else if (name == "ConversionMassSpecType") {
            if (options.find("-instrument") == options.end()) {
                options["-instrument"] = parameter;
            }
        } else if (name == "ConversionWarpedMesh" && parameter == "1") {
            if (options.find("-warped") == options.end()) {
                options["-warped"] = parameter;
            }
        } else {
            // ignoring unknown parameters...
        }
    }
    return false;
}

void print_parameters_summary(const Grid::Parameters& parameters) {
    std::cout << "The following parameters were set:" << std::endl;
    // Dimensions.
    std::cout << "DIMENSIONS:" << std::endl;
    std::cout << "num_mz:" << parameters.dimensions.n << std::endl;
    std::cout << "num_rt:" << parameters.dimensions.m << std::endl;
    // Bounds.
    std::cout << "BOUNDS:" << std::endl;
    std::cout << "min_rt:" << parameters.bounds.min_rt << std::endl;
    std::cout << "max_rt:" << parameters.bounds.max_rt << std::endl;
    std::cout << "min_mz:" << parameters.bounds.min_mz << std::endl;
    std::cout << "max_mz:" << parameters.bounds.max_mz << std::endl;
    // SmoothingParams.
    std::cout << "SMOOTHING PARAMETERS:" << std::endl;
    std::cout << "mz:" << parameters.smoothing_params.mz << std::endl;
    std::cout << "sigma_mz:" << parameters.smoothing_params.sigma_mz
              << std::endl;
    std::cout << "sigma_rt:" << parameters.smoothing_params.sigma_rt
              << std::endl;
    // Instrument type.
    std::cout << "INSTRUMENT TYPE:" << std::endl;
    std::cout << parameters.instrument_type << std::endl;

    // Flags.
    std::cout << "FLAGS:" << std::endl;
    std::cout << "Warped: " << bool(parameters.flags & Grid::Flags::WARPED_MESH)
              << std::endl;

    // Memory usage.
    double x = parameters.dimensions.n * parameters.dimensions.m * 8;
    std::cout << "APPROXIMATE MEMORY USAGE (BYTES):" << x << std::endl;
}

int main(int argc, char* argv[]) {
    // Flag format is map where the key is the flag name and contains a tuple
    // with the description and if it takes extra parameters or not:
    // <description, takes_parameters>
    const std::map<std::string, std::pair<std::string, bool>> accepted_flags = {
        // Dimensions.
        {"-num_mz", {"The number of sampling points for the grid on mz", true}},
        {"-num_rt", {"The number of sampling points for the grid on rt", true}},
        {"-delta_mz",
         {"The interval between sampling points for the grid on mz", true}},
        {"-delta_rt",
         {"The interval between sampling points for the grid on mz", true}},
        // Bounds.
        {"-min_rt", {"The minimum rt value", true}},
        {"-max_rt", {"The maximum rt value", true}},
        {"-min_mz", {"The minimum mz value", true}},
        {"-max_mz", {"The maximum mz value", true}},
        // SmoothingParams.
        {"-smooth_mz",
         {"The mass at which the smoothing sigma is given", true}},
        {"-sigma_mz", {"The smoothing sigma in the mz direction", true}},
        {"-sigma_rt", {"The smoothing sigma in the rt direction", true}},
        // Instrument::Type.
        {"-instrument",
         {"The instrument in which the data was extracted", true}},
        // Flags.
        {"-warped", {"Specify if the output grid will be warped", false}},
        // Command parameters.
        {"-out_dir", {"The output directory", true}},
        {"-help", {"Display available options", false}},
        {"-config", {"Specify the configuration file", true}},
    };

    if (argc == 1) {
        print_usage();
        return -1;
    }

    // Parse arguments and extract options and file names.
    options_map options;
    std::vector<std::string> files;
    for (size_t i = 1; i < argc; ++i) {
        auto arg = argv[i];
        if (arg[0] == '-') {
            if (accepted_flags.find(arg) == accepted_flags.end()) {
                std::cout << "unknown option: " << arg << std::endl;
                print_usage();
                return -1;
            }

            auto flag = accepted_flags.at(arg);
            if (flag.second) {
                if (i + 1 >= argc || argv[i + 1][0] == '-') {
                    std::cout << "no parameters specified for " << arg
                              << std::endl;
                    print_usage();
                    return -1;
                }
                ++i;
                options[arg] = argv[i];
            } else {
                options[arg] = "";
            }
        } else {
            files.emplace_back(arg);
        }
    }

    if (options.find("-help") != options.end()) {
        print_usage();
        // Find maximum option length to adjust text padding.
        auto padding = 0;
        for (const auto& e : accepted_flags) {
            if (e.first.size() > padding) {
                padding = e.first.size();
            }
        }

        // Print options with a 4 space padding between flag name and
        // description.
        std::cout << "OPTIONS:" << std::endl;
        for (const auto& e : accepted_flags) {
            std::cout << e.first;
            // If the option requires an argument we have to specify it,
            // otherwise we add padding.
            if (e.second.second) {
                std::cout << " <arg>";
            } else {
                std::cout << "      ";
            }
            for (size_t i = 0; i < (padding - e.first.size()) + 4; ++i) {
                std::cout << " ";
            }
            std::cout << e.second.first << std::endl;
        }
        return 0;
    }

    if (files.empty()) {
        std::cout << "No input files specified." << std::endl;
        print_usage();
        return -1;
    }

    // If config file is provided, read it and parse it. The parameters
    // specified as command line arguments will override the config file.
    if (options.find("-config") != options.end()) {
        std::filesystem::path config_path = options["-config"];
        if (!std::filesystem::exists(config_path)) {
            std::cout << "error: couldn't find config file " << config_path
                      << std::endl;
            print_usage();
            return -1;
        }

        // TODO(alex): accept both json and hdr for now.
        auto extension = config_path.extension();
        if (extension == ".json") {
            // TODO(alex): parse .json file...
        } else if (extension == ".hdr") {
            parse_hdr(config_path, options);
        } else {
            std::cout << "error: invalid format for config file " << config_path
                      << std::endl;
            print_usage();
            return -1;
        }
    }

    // Parse the options to build the Grid::Parameters struct.
    Grid::Parameters parameters;

    // Get the bounds.
    if ((options.find("-min_rt") == options.end()) ||
        (options.find("-max_rt") == options.end()) ||
        (options.find("-min_mz") == options.end()) ||
        (options.find("-max_mz") == options.end())) {
        std::cout
            << "Grid bounds (min_rt, max_rt, min_mz, max_mz) not specified"
            << std::endl;
        return -1;
    }
    auto min_rt = options["-min_rt"];
    auto max_rt = options["-max_rt"];
    auto min_mz = options["-min_mz"];
    auto max_mz = options["-max_mz"];
    if (!is_number(min_rt)) {
        std::cout << "error: "
                  << "min_rt"
                  << " has to be a number" << std::endl;
        print_usage();
        return -1;
    }
    if (!is_number(max_rt)) {
        std::cout << "error: "
                  << "max_rt"
                  << " has to be a number" << std::endl;
        print_usage();
        return -1;
    }
    if (!is_number(min_mz)) {
        std::cout << "error: "
                  << "min_mz"
                  << " has to be a number" << std::endl;
        print_usage();
        return -1;
    }
    if (!is_number(max_mz)) {
        std::cout << "error: "
                  << "max_mz"
                  << " has to be a number" << std::endl;
        print_usage();
        return -1;
    }
    parameters.bounds.min_rt = std::stod(min_rt);
    parameters.bounds.max_rt = std::stod(max_rt);
    parameters.bounds.min_mz = std::stod(min_mz);
    parameters.bounds.max_mz = std::stod(max_mz);

    // Get the dimensions. The number of sampling points in either direction can
    // be set manually specifying -num_mz and -num_rt or indirectly using
    // -delta_mz and -delta_rt. The former will take priority over the later.
    if ((options.find("-delta_mz") != options.end()) &&
        (options.find("-num_mz") == options.end())) {
        auto delta_mz = options["-delta_mz"];
        if (!is_number(delta_mz)) {
            std::cout << "error: "
                      << "delta_mz"
                      << " has to be a number" << std::endl;
            print_usage();
            return -1;
        }
        unsigned int num_mz =
            (parameters.bounds.max_mz - parameters.bounds.min_mz) /
            std::stod(delta_mz);
        options["-num_mz"] = std::to_string(num_mz);
    }
    if ((options.find("-delta_rt") != options.end()) &&
        (options.find("-num_rt") == options.end())) {
        auto delta_rt = options["-delta_rt"];
        if (!is_number(delta_rt)) {
            std::cout << "error: "
                      << "delta_rt"
                      << " has to be a number" << std::endl;
            print_usage();
            return -1;
        }
        unsigned int num_rt =
            (parameters.bounds.max_rt - parameters.bounds.min_rt) /
            std::stod(delta_rt);
        options["-num_rt"] = std::to_string(num_rt);
    }
    if ((options.find("-num_mz") == options.end()) ||
        (options.find("-num_rt") == options.end())) {
        std::cout << "Grid dimensions (num_mz, num_rt) not specified"
                  << std::endl;
        return -1;
    }
    auto num_mz = options["-num_mz"];
    auto num_rt = options["-num_rt"];
    if (!is_unsigned_int(num_mz)) {
        std::cout << "error: "
                  << "num_mz"
                  << " has to be a positive integer" << std::endl;
        print_usage();
        return -1;
    }
    if (!is_unsigned_int(num_rt)) {
        std::cout << "error: "
                  << "num_rt"
                  << " has to be a positive integer" << std::endl;
        print_usage();
        return -1;
    }
    parameters.dimensions.n = std::stoi(num_mz);
    parameters.dimensions.m = std::stoi(num_rt);

    // Get the smoothing parameters.
    if ((options.find("-smooth_mz") == options.end()) ||
        (options.find("-sigma_mz") == options.end()) ||
        (options.find("-sigma_rt") == options.end())) {
        std::cout << "Smoothing parameters (smooth_mz, sigma_mz, sigma_rt) not "
                     "specified"
                  << std::endl;
        return -1;
    }
    auto smooth_mz = options["-smooth_mz"];
    auto sigma_mz = options["-sigma_mz"];
    auto sigma_rt = options["-sigma_rt"];
    if (!is_number(smooth_mz)) {
        std::cout << "error: "
                  << "smooth_mz"
                  << " has to be a number" << std::endl;
        print_usage();
        return -1;
    }
    if (!is_number(sigma_mz)) {
        std::cout << "error: "
                  << "sigma_mz"
                  << " has to be a number" << std::endl;
        print_usage();
        return -1;
    }
    if (!is_number(sigma_rt)) {
        std::cout << "error: "
                  << "sigma_rt"
                  << " has to be a number" << std::endl;
        print_usage();
        return -1;
    }
    parameters.smoothing_params.mz = std::stod(smooth_mz);
    parameters.smoothing_params.sigma_mz = std::stod(sigma_mz);
    parameters.smoothing_params.sigma_rt = std::stod(sigma_rt);

    // Get the instrument type.
    if (options.find("-instrument") == options.end()) {
        std::cout << "Instrument type (instrument) not specified" << std::endl;
        print_usage();
        return -1;
    }
    auto instrument = options["-instrument"];
    // Transform instrument to lowercase to prevent typos.
    for (auto& ch : instrument) {
        ch = std::tolower(ch);
    }
    if (instrument == "orbitrap") {
        parameters.instrument_type = Instrument::Type::ORBITRAP;
    } else if (instrument == "quad" || instrument == "iontrap") {
        parameters.instrument_type = Instrument::Type::QUAD;
    } else if (instrument == "tof" || instrument == "qtof") {
        parameters.instrument_type = Instrument::Type::TOF;
    } else if (instrument == "fticr") {
        parameters.instrument_type = Instrument::Type::FTICR;
    } else {
        std::cout << "Unknown instrument type: " << instrument << std::endl;
        return -1;
    }

    // Set up the flags.
    if (options.find("-warped") != options.end()) {
        parameters.flags |= Grid::Flags::WARPED_MESH;
    }

    // Set up the output directory and check if it exists.
    if (options.find("-out_dir") == options.end()) {
        options["-out_dir"] = ".";
    }
    if (!std::filesystem::exists(options["-out_dir"])) {
        std::cout << "error: couldn't find output directory \""
                  << options["-out_dir"] << "\"" << std::endl;
        print_usage();
        return -1;
    }

    // Execute the program here.
    for (const auto& file_name : files) {
        std::filesystem::path input_file = file_name;
        // Check if the files exist.
        if (!std::filesystem::exists(input_file)) {
            std::cout << "error: couldn't find file " << input_file << std::endl;
            print_usage();
            return -1;
        }

        // Check if the file has the appropriate format.
        std::string extension = input_file.extension();
        std::string lowercase_extension = extension;
        for (auto& ch : lowercase_extension) {
            ch = std::tolower(ch);
        }
        if (lowercase_extension == ".mzxml") {
            print_parameters_summary(parameters);

            // Instantiate memory.
            std::vector<double> data(parameters.dimensions.n *
                                     parameters.dimensions.m);
            // Open input file.
            std::ifstream stream;
            stream.open(input_file);
            if (!stream) {
                std::cout << "error: could not open input file " << input_file
                          << std::endl;
                return -1;
            }

            // Prepare the name of the output file.
            auto datfile_name = options["-out_dir"] /
                                input_file.filename().replace_extension(".dat");
            std::ofstream datfile;
            datfile.open(datfile_name, std::ios::out | std::ios::binary);
            if (!datfile) {
                std::cout << "error: could not open file " << datfile_name
                          << " for writing" << std::endl;
                return -1;
            }

            std::cout << "Parsing file..." << std::endl;
            auto scans = XmlReader::read_next_scan(stream, parameters);
            if (scans == std::nullopt) {
                std::cout << "error: no scans found on file " << input_file
                          << " for the given parameters" << std::endl;
                return -1;
            }
            do {
                for (const auto& scan : scans.value()) {
                    Grid::splat(scan.mz, scan.rt, scan.value, parameters, data);
                }
                scans = XmlReader::read_next_scan(stream, parameters);
            } while (scans != std::nullopt);

            std::cout << "Saving into dat file..." << std::endl;
            Grid::File::write(datfile, data, parameters);
        } else {
            std::cout << "error: unknown file format for file " << input_file
                      << std::endl;
            print_usage();
            return -1;
        }
    }

    return 0;
}
