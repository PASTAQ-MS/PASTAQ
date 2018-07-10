#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <thread>
#include <vector>

#include "grid.hpp"
#include "grid_files.hpp"
#include "grid_runners.hpp"
#include "xml_reader.hpp"

// Type aliases.
using options_map = std::map<std::string, std::string>;

void print_usage() {
    std::cout << "USAGE: grid [-help] [options] <files>" << std::endl;
}

// Helper functions to check if the given string contains a number.
bool is_unsigned_int(std::string &s) {
    std::regex double_regex("^([[:digit:]]+)$");
    return std::regex_search(s, double_regex);
}
bool is_number(std::string &s) {
    std::regex double_regex("^([[:digit:]]+[\\.]?[[:digit:]]*)$");
    return std::regex_search(s, double_regex);
}

// Helper function to trim the whitespace surrounding a string.
void trim_space(std::string &s) {
    auto not_space = [](int ch) { return !std::isspace(ch); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
    s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
}

void parse_json(const std::filesystem::path &path, options_map &options,
                std::vector<std::string> &files) {
    std::ifstream stream(path);
    std::string line;
    std::string content;

    // Read all the lines.
    while (stream.good()) {
        std::getline(stream, line);
        trim_space(line);
        if (line.empty() || line[0] == '#') {
            continue;
        }
        content += line + " ";
    }

    // Find the contents of the grid configuration.
    std::regex grid_regex(R"("grid"[[:space:]]*:[[:space:]]*\{(.*)\})");
    std::smatch matches;
    std::regex_search(content, matches, grid_regex);
    if (matches.size() != 2 || matches[1].str().empty()) {
        std::cout << "error: could not find \"grid\" on the config file"
                  << std::endl;
        std::exit(-1);
    }
    content = matches[1];

    // Match the different parts of the grid configuration. Note that we are
    // performing a very simplistic parsing, the right thing to do would be to
    // create an AST and walk it to perform the parsing. This simplistic
    // approach means that if the user mistakenly puts the "paths" inside
    // "parameters" or "config" it will still be accepted properly.
    //
    // Furthermore, we are not performing any validation on the JSON file, so a
    // malformed file in principle could be accepted.
    {
        // Parse file paths.
        auto pos = content.find("\"paths\"");
        if (pos != std::string::npos) {
            auto begin = content.find('[', pos) + 1;
            auto end = content.find(']', begin) - 1;
            if (begin == std::string::npos || end == std::string::npos ||
                end < begin) {
                std::cout << "error: malformed config file" << std::endl;
                std::exit(-1);
            }

            auto config_files = content.substr(begin, end - begin);
            for (auto &ch : config_files) {
                if (ch == ',' || ch == '"') {
                    ch = ' ';
                }
            }
            std::stringstream ss(config_files);
            while (ss.good()) {
                std::string file_name;
                ss >> file_name;
                if (!file_name.empty() && std::find(files.begin(), files.end(),
                                                    file_name) == files.end()) {
                    files.push_back(file_name);
                }
            }
        }
    }

    {
        // Parse parameters.
        auto pos = content.find("\"parameters\"");
        if (pos != std::string::npos) {
            auto begin = content.find('{', pos) + 1;
            auto end = content.find('}', begin) - 1;
            if (begin == std::string::npos || end == std::string::npos ||
                end < begin) {
                std::cout << "error: malformed config file" << std::endl;
                std::exit(-1);
            }

            auto config_parameters = content.substr(begin, end - begin);
            for (auto &ch : config_parameters) {
                if (ch == ',' || ch == ':' || ch == '"') {
                    ch = ' ';
                }
            }

            std::stringstream ss(config_parameters);
            while (ss.good()) {
                std::string name;
                std::string value;
                ss >> name;
                ss >> value;
                name.insert(0, "-");
                if (options.find(name) == options.end()) {
                    options[name] = value;
                }
            }
        }
    }

    {
        // Parse config.
        auto pos = content.find("\"config\"");
        if (pos != std::string::npos) {
            auto begin = content.find('{', pos) + 1;
            auto end = content.find('}', begin) - 1;
            if (begin == std::string::npos || end == std::string::npos ||
                end < begin) {
                std::cout << "error: malformed config file" << std::endl;
                std::exit(-1);
            }

            auto config_config = content.substr(begin, end - begin);
            for (auto &ch : config_config) {
                if (ch == ',' || ch == ':' || ch == '"') {
                    ch = ' ';
                }
            }

            std::stringstream ss(config_config);
            while (ss.good()) {
                std::string name;
                std::string value;
                ss >> name;
                ss >> value;
                name.insert(0, "-");
                if (options.find(name) == options.end()) {
                    options[name] = value;
                }
            }
        }
    }
}

bool parse_hdr(const std::filesystem::path &path, options_map &options) {
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
        } else if (name == "ConversionMassSpecType") {
            if (options.find("-instrument") == options.end()) {
                options["-instrument"] = parameter;
            }
        } else if (name == "ConversionWarpedMesh" && parameter == "1") {
            if (options.find("-warped") == options.end()) {
                options["-warped"] = "true";
            }
        } else {
            // ignoring unknown parameters...
        }
    }
    return false;
}

void print_parameters_summary(const Grid::Parameters &parameters) {
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

int main(int argc, char *argv[]) {
    // Flag format is map where the key is the flag name and contains a tuple
    // with the description and if it takes extra parameters or not:
    // <description, takes_parameters>
    const std::map<std::string, std::pair<std::string, bool>> accepted_flags = {
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
        {"-rawdump",
         {"Enable the dump of the raw points inside the given bounds", false}},
        {"-parallel", {"Enable parallel processing", false}},
        {"-n_threads",
         {"Specify the maximum number of threads that will be used for the "
          "calculations",
          true}},
    };

    if (argc == 1) {
        print_usage();
        return -1;
    }

    // Parse arguments and extract options and file names.
    options_map options;
    std::vector<std::string> files;
    for (int i = 1; i < argc; ++i) {
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
        size_t padding = 0;
        for (const auto &e : accepted_flags) {
            if (e.first.size() > padding) {
                padding = e.first.size();
            }
        }

        // Print options with a 4 space padding between flag name and
        // description.
        std::cout << "OPTIONS:" << std::endl;
        for (const auto &e : accepted_flags) {
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

        // FIXME(alex): accept both json and hdr for now. Ideally we would stick
        // with just one configuration format.
        auto extension = config_path.extension();
        if (extension == ".json") {
            parse_json(config_path, options, files);
        } else if (extension == ".hdr") {
            parse_hdr(config_path, options);
        } else {
            std::cout << "error: invalid format for config file " << config_path
                      << std::endl;
            print_usage();
            return -1;
        }
    }

    if (files.empty()) {
        std::cout << "No input files specified." << std::endl;
        print_usage();
        return -1;
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
    for (auto &ch : instrument) {
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

    if (options.find("-warped") != options.end() &&
        (options["-warped"] == "true" || options["-warped"].empty())) {
        // Set up the flags.
        parameters.flags |= Grid::Flags::WARPED_MESH;
    }

    // To calculate the dimensions of the grid we need to use the bounds and
    // the reference sigma_mz at a given mass for the given instrument.
    Grid::calculate_dimensions(parameters);

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

    // Set up maximum concurrency.
    unsigned int max_threads = std::thread::hardware_concurrency();
    if (!max_threads) {
        std::cout << "error: this system does not support parallel processing"
                  << std::endl;
        return -1;
    }
    if ((options.find("-n_threads") != options.end()) &&
        ((options.find("-parallel") != options.end()) &&
         (options["-parallel"] == "true" || options["-parallel"].empty()))) {
        auto n_threads = options["-n_threads"];
        if (!is_unsigned_int(n_threads)) {
            std::cout << "error: "
                      << "n_threads"
                      << " has to be a positive integer" << std::endl;
            print_usage();
            return -1;
        }
        max_threads = std::stoi(n_threads);
    }

    // Execute the program here.
    for (const auto &file_name : files) {
        std::filesystem::path input_file = file_name;
        // Check if the files exist.
        if (!std::filesystem::exists(input_file)) {
            std::cout << "error: couldn't find file " << input_file
                      << std::endl;
            print_usage();
            return -1;
        }

        // Check if the file has the appropriate format.
        std::string extension = input_file.extension();
        std::string lowercase_extension = extension;
        for (auto &ch : lowercase_extension) {
            ch = std::tolower(ch);
        }
        if (lowercase_extension == ".mzxml") {
            print_parameters_summary(parameters);

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
            std::ofstream datfile_stream;
            datfile_stream.open(datfile_name, std::ios::out | std::ios::binary);
            if (!datfile_stream) {
                std::cout << "error: could not open file " << datfile_name
                          << " for writing" << std::endl;
                return -1;
            }

            std::cout << "Parsing file..." << std::endl;
            auto points = XmlReader::read_next_scan(stream, parameters);
            if (points == std::nullopt) {
                std::cout << "error: no points found on file " << input_file
                          << " for the given parameters" << std::endl;
                return -1;
            }
            std::vector<Grid::Point> all_points = {};
            do {
                all_points.insert(end(all_points), begin(points.value()),
                                  end(points.value()));
                points = XmlReader::read_next_scan(stream, parameters);
            } while (points != std::nullopt);

            if (options.find("-rawdump") != options.end() &&
                (options["-rawdump"] == "true" ||
                 options["-rawdump"].empty())) {
                std::cout << "Generating rawdump file..." << std::endl;
                // Prepare the name of the rawdump output file.
                auto rawdump_name =
                    options["-out_dir"] /
                    input_file.filename().replace_extension(".rawdump");
                std::ofstream rawdump_stream;
                rawdump_stream.open(rawdump_name,
                                    std::ios::out | std::ios::binary);
                if (!rawdump_stream) {
                    std::cout << "error: could not open file " << rawdump_name
                              << " for writing" << std::endl;
                    return -1;
                }

                // Write rawdump to disk.
                if (!Grid::Files::Rawdump::write(rawdump_stream, all_points)) {
                    std::cout
                        << "error: the raw dump could not be saved properly"
                        << std::endl;
                    return -1;
                }
            }

            std::cout << "Performing grid splatting..." << std::endl;
            std::vector<double> data;
            if ((options.find("-parallel") != options.end()) &&
                (options["-parallel"] == "true" ||
                 options["-parallel"].empty())) {
                data = Grid::Runners::Parallel::run(max_threads, parameters,
                                                    all_points);
            } else {
                data = Grid::Runners::Serial::run(parameters, all_points);
            }

            std::cout << "Saving grid into dat file..." << std::endl;
            if (!Grid::Files::Dat::write(datfile_stream, data, parameters)) {
                std::cout << "error: the grid could not be saved properly"
                          << std::endl;
                return -1;
            }
        } else if (lowercase_extension == ".rawdump") {
            print_parameters_summary(parameters);

            // Prepare the name of the output file.
            auto datfile_name = options["-out_dir"] /
                                input_file.filename().replace_extension(".dat");
            std::ofstream datfile_stream;
            datfile_stream.open(datfile_name, std::ios::out | std::ios::binary);
            if (!datfile_stream) {
                std::cout << "error: could not open file " << datfile_name
                          << " for writing" << std::endl;
                return -1;
            }

            // Open the file for reading.
            std::ifstream stream;
            stream.open(input_file, std::ios::in | std::ios::binary);
            if (!stream) {
                std::cout << "error: could not open input file " << input_file
                          << std::endl;
                return -1;
            }

            // Load the points into memory.
            std::vector<Grid::Point> all_points;
            if (!Grid::Files::Rawdump::read(stream, all_points)) {
                std::cout << "error: the raw dump could not be loaded properly"
                          << std::endl;
                return -1;
            }
            if (all_points.empty()) {
                std::cout << "error: the raw dump does not contain any points"
                          << std::endl;
                return -1;
            }
            std::cout << "Loaded " << all_points.size() << " points"
                      << std::endl;

            std::cout << "Performing grid splatting..." << std::endl;
            std::vector<double> data;
            if ((options.find("-parallel") != options.end()) &&
                (options["-parallel"] == "true" ||
                 options["-parallel"].empty())) {
                data = Grid::Runners::Parallel::run(max_threads, parameters,
                                                    all_points);
            } else {
                data = Grid::Runners::Serial::run(parameters, all_points);
            }

            std::cout << "Saving grid into dat file..." << std::endl;
            if (!Grid::Files::Dat::write(datfile_stream, data, parameters)) {
                std::cout << "error: the grid could not be saved properly"
                          << std::endl;
                return -1;
            }
        } else {
            std::cout << "error: unknown file format for file " << input_file
                      << std::endl;
            print_usage();
            return -1;
        }
    }

    return 0;
}