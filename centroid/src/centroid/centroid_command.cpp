#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <thread>
#include <vector>

#include "centroid/centroid.hpp"
#include "centroid/centroid_files.hpp"
#include "centroid/centroid_runners.hpp"
#include "grid/grid_files.hpp"

// Type aliases.
using options_map = std::map<std::string, std::string>;

void print_usage() {
    std::cout << "USAGE: centroid [-help] [options] <files>" << std::endl;
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

    // Find the contents of the centroid configuration.
    std::regex centroid_regex(R"("centroid"[[:space:]]*:[[:space:]]*\{(.*)\})");
    std::smatch matches;
    std::regex_search(content, matches, centroid_regex);
    if (matches.size() != 2 || matches[1].str().empty()) {
        std::cout << "error: could not find \"centroid\" on the config file"
                  << std::endl;
        std::exit(-1);
    }
    content = matches[1];

    // Match the different parts of the centroid configuration. Note that we are
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

        if (name == "ConversionPeakThreshold") {
            if (options.find("-threshold") == options.end()) {
                options["-threshold"] = parameter;
            }
        } else if (name == "ConversionNPeaksToFind") {
            if (options.find("-n_peaks") == options.end()) {
                options["-n_peaks"] = parameter;
            }
        } else {
            // ignoring unknown parameters...
        }
    }
    return false;
}

int main(int argc, char *argv[]) {
    // Flag format is map where the key is the flag name and contains a tuple
    // with the description and if it takes extra parameters or not:
    // <description, takes_parameters>
    const std::map<std::string, std::pair<std::string, bool>> accepted_flags = {
        // Centroid parameters.
        {"-threshold",
         {"The minimum height value considered to be part of a peak", true}},
        {"-n_peaks",
         {"The maximum number of peaks to be reported. Filter peaks based on "
          "descending local max height. (If set to 0, no peaks will be "
          "filtered)",
          true}},
        // Command parameters.
        {"-out_dir", {"The output directory", true}},
        {"-help", {"Display available options", false}},
        {"-config", {"Specify the configuration file", true}},
        {"-csvdump",
         {"Dump the detected peaks as a `csv` file in addition to the `.bpks`",
          false}},
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

    // Parse the options to build the Centroid::Parameters struct.
    Centroid::Parameters parameters = {};
    if (options.find("-n_peaks") != options.end()) {
        auto n_peaks = options["-n_peaks"];
        if (!is_unsigned_int(n_peaks)) {
            std::cout << "error: n_peaks has to be a positive integer"
                      << std::endl;
            print_usage();
            return -1;
        }
        parameters.n_peaks = std::stoi(n_peaks);
    }
    if (options.find("-threshold") != options.end()) {
        auto threshold = options["-threshold"];
        if (!is_number(threshold)) {
            std::cout << "error: threshold has to be a positive integer"
                      << std::endl;
            print_usage();
            return -1;
        }
        parameters.threshold = std::stod(threshold);
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

    // Set up maximum concurrency.
    uint64_t max_threads = std::thread::hardware_concurrency();
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
            std::cout << "error: n_threads has to be a positive integer"
                      << std::endl;
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
        if (lowercase_extension == ".dat") {
            // Open input file.
            std::ifstream stream;
            stream.open(input_file);
            if (!stream) {
                std::cout << "error: could not open input file " << input_file
                          << std::endl;
                return -1;
            }

            // Prepare the name of the output file.
            auto bpks_name = options["-out_dir"] /
                             input_file.filename().replace_extension(".bpks");
            std::ofstream bpks_stream;
            bpks_stream.open(bpks_name, std::ios::out | std::ios::binary);
            if (!bpks_stream) {
                std::cout << "error: could not open file " << bpks_name
                          << " for writing" << std::endl;
                return -1;
            }

            std::cout << "Loading file..." << std::endl;
            std::vector<double> grid_data;
            if (!Grid::Files::Dat::read(stream, &grid_data,
                                        &parameters.grid_params)) {
                std::cout << "error: loading the data from " << input_file
                          << std::endl;
                return -1;
            }

            std::vector<Centroid::Peak> peaks;
            std::cout << "Building peaks..." << std::endl;
            if ((options.find("-parallel") != options.end()) &&
                (options["-parallel"] == "true" ||
                 options["-parallel"].empty())) {
                peaks = Centroid::Runners::Parallel::run(max_threads,
                                                         parameters, grid_data);
            } else {
                peaks = Centroid::Runners::Serial::run(parameters, grid_data);
            }
            if (peaks.size() == 0) {
                std::cout << "error: no peaks could be found on " << input_file
                          << std::endl;
                return -1;
            }

            // TODO(alex): Is this necessary? We are already sorting the peaks
            // by the detected local max.
            std::cout << "Sorting peaks by height (centroid)..." << std::endl;
            auto sort_peaks = [](const Centroid::Peak &p1,
                                 const Centroid::Peak &p2) -> bool {
                return (p1.height_centroid > p2.height_centroid) ||
                       ((p1.height_centroid == p2.height_centroid) &&
                        (p1.total_intensity_centroid >
                         p2.total_intensity_centroid));
            };
            std::stable_sort(peaks.begin(), peaks.end(), sort_peaks);

            std::cout << "Saving peaks to bpks file..." << std::endl;
            if (!Centroid::Files::Bpks::write_peaks(
                    bpks_stream, parameters.grid_params, peaks)) {
                std::cout << "error: the peaks could not be saved properly"
                          << std::endl;
                return -1;
            }

            if (options.find("-csvdump") != options.end() &&
                (options["-csvdump"] == "true" ||
                 options["-csvdump"].empty())) {
                std::cout << "Dumping peaks to csv..." << std::endl;
                // Prepare the name of the csvdump output file.
                auto csvdump_name =
                    options["-out_dir"] /
                    input_file.filename().replace_extension(".csv");
                std::ofstream csvdump_stream;
                csvdump_stream.open(csvdump_name,
                                    std::ios::out | std::ios::binary);
                if (!csvdump_stream) {
                    std::cout << "error: could not open file " << csvdump_name
                              << " for writing" << std::endl;
                    return -1;
                }

                // Write csv file to disk.
                if (!Centroid::Files::Csv::write_peaks(csvdump_stream, peaks)) {
                    std::cout
                        << "error: the csv dump could not be saved properly"
                        << std::endl;
                    return -1;
                }
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
