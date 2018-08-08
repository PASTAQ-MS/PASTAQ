#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <thread>
#include <vector>

#include "centroid/centroid_files.hpp"
#include "grid/grid.hpp"
#include "warp2d/warp2d.hpp"
#include "warp2d/warp2d_runners.hpp"

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
    std::regex grid_regex(R"("warp2d"[[:space:]]*:[[:space:]]*\{(.*)\})");
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

int main(int argc, char *argv[]) {
    // Flag format is map where the key is the flag name and contains a tuple
    // with the description and if it takes extra parameters or not:
    // <description, takes_parameters>
    const std::map<std::string, std::pair<std::string, bool>> accepted_flags = {
        // Warp2D parameters.
        {"-reference_file",
         {"The path to the file to be used as a reference", true}},
        {"-slack",
         {"The maximum number of time points that can be warped", true}},
        {"-window_size",
         {"The size of the window that will be used for segmenting the data",
          true}},
        {"-num_points",
         {"The number of time points in which the retention time range will be "
          "divided",
          true}},
        // Command parameters.
        {"-out_dir", {"The output directory", true}},
        {"-help", {"Display available options", false}},
        {"-config", {"Specify the configuration file", true}},
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
    Warp2D::Parameters parameters = {};

    // Get the Warp2D parameters.
    if (options.find("-slack") == options.end()) {
        std::cout << "Slack (slack) not specified" << std::endl;
        return -1;
    }
    auto slack = options["-slack"];
    if (!is_unsigned_int(slack)) {
        std::cout << "error: slack has to be a positive integer" << std::endl;
        print_usage();
        return -1;
    }
    parameters.slack = std::stoi(slack);

    if (options.find("-window_size") == options.end()) {
        std::cout << "Window size (window_size) not specified" << std::endl;
        return -1;
    }
    auto window_size = options["-window_size"];
    if (!is_unsigned_int(window_size)) {
        std::cout << "error: window_size has to be a positive integer"
                  << std::endl;
        print_usage();
        return -1;
    }
    parameters.window_size = std::stoi(window_size);

    if (options.find("-num_points") == options.end()) {
        std::cout << "Number of rt points (num_points) not specified"
                  << std::endl;
        return -1;
    }
    auto num_points = options["-num_points"];
    if (!is_unsigned_int(num_points)) {
        std::cout << "error: num_points has to be a positive integer"
                  << std::endl;
        print_usage();
        return -1;
    }
    parameters.num_points = std::stoi(num_points);

    if (options.find("-reference_file") == options.end()) {
        std::cout << "Reference file (reference_file) not specified"
                  << std::endl;
        return -1;
    }
    auto reference_file = options["-reference_file"];

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

    // Read the reference file peaks.
    std::cout << "Reading the reference file peaks" << std::endl;
    std::vector<Centroid::Peak> reference_peaks;
    Grid::Parameters reference_grid_params;
    {
        std::ifstream stream;
        stream.open(reference_file);
        if (!stream) {
            std::cout << "error: could not open reference file "
                      << reference_file << std::endl;
            return -1;
        }
        std::filesystem::path input_file = reference_file;
        // Check if the file has the appropriate format.
        std::string extension = input_file.extension();
        std::string lowercase_extension = extension;
        for (auto &ch : lowercase_extension) {
            ch = std::tolower(ch);
        }

        if (lowercase_extension == ".bpks") {
            // Read into peak array.
            if (!Centroid::Files::Bpks::read_peaks(
                    stream, &reference_grid_params, &reference_peaks)) {
                std::cout
                    << "error: couldn't read peaks from the reference file"
                    << std::endl;
                return -1;
            }
        } else if (lowercase_extension == ".csv") {
            // TODO: read into peak array.
            // Centroid::Files::Csv::read_peaks(stream, &grid_params,
            // &peaks[0]);
        } else {
            std::cout << "error: unknown file format for file " << input_file
                      << std::endl;
            print_usage();
            return -1;
        }
    }

    // Execute the program here.
    for (const auto &file_name : files) {
        std::filesystem::path input_file = file_name;

        // Open input file.
        std::ifstream stream;
        stream.open(input_file);
        if (!stream) {
            std::cout << "error: could not open input file " << input_file
                      << std::endl;
            return -1;
        }

        // Check if the file has the appropriate format.
        std::string extension = input_file.extension();
        std::string lowercase_extension = extension;
        for (auto &ch : lowercase_extension) {
            ch = std::tolower(ch);
        }

        std::vector<Centroid::Peak> peaks;
        Grid::Parameters grid_params;
        if (lowercase_extension == ".bpks") {
            // Read into peak array.
            std::cout << "Reading peaks from file: " << input_file << std::endl;
            if (!Centroid::Files::Bpks::read_peaks(stream, &grid_params,
                                                   &peaks)) {
                std::cout
                    << "error: couldn't read peaks from the reference file"
                    << std::endl;
                return -1;
            }
        } else if (lowercase_extension == ".csv") {
            // TODO: read into peak array.
            // Centroid::Files::Csv::read_peaks(stream, &grid_params,
            // &peaks[0]);
        } else {
            std::cout << "error: unknown file format for file " << input_file
                      << std::endl;
            print_usage();
            return -1;
        }

        // Prepare the name of the output file.
        auto outfile_name = options["-out_dir"] / input_file.filename();
        std::ofstream outfile_stream;
        outfile_stream.open(outfile_name, std::ios::out | std::ios::binary);
        if (!outfile_stream) {
            std::cout << "error: could not open file " << outfile_name
                      << " for writing" << std::endl;
            return -1;
        }
        // Perform warping.
        std::cout << "Performing warping..." << std::endl;
        std::vector<Centroid::Peak> warped_peaks;
        if ((options.find("-parallel") != options.end()) &&
            (options["-parallel"] == "true" || options["-parallel"].empty())) {
            warped_peaks = Warp2D::Runners::Parallel::run(
                reference_peaks, peaks, parameters, max_threads);
        } else {
            warped_peaks = Warp2D::Runners::Serial::run(reference_peaks, peaks,
                                                        parameters);
        }

        if (lowercase_extension == ".bpks") {
            std::cout << "Saving peaks to disk in bpks..." << std::endl;
            if (!Centroid::Files::Bpks::write_peaks(outfile_stream, grid_params,
                                                    warped_peaks)) {
                std::cout << "error: couldn't write warped peaks into file "
                          << outfile_name << std::endl;
                return -1;
            }
            // DEBUG
            std::cout << "Saving peaks to disk in csv..." << std::endl;
            auto csv_outfile_name =
                options["-out_dir"] /
                input_file.filename().replace_extension(".csv");
            std::ofstream csv_outfile_stream;
            csv_outfile_stream.open(csv_outfile_name,
                                    std::ios::out | std::ios::binary);
            if (!Centroid::Files::Csv::write_peaks(csv_outfile_stream,
                                                   warped_peaks)) {
                std::cout << "error: couldn't write warped peaks into file "
                          << csv_outfile_name << std::endl;
                return -1;
            }
        } else if (lowercase_extension == ".csv") {
            if (!Centroid::Files::Csv::write_peaks(outfile_stream,
                                                   warped_peaks)) {
                std::cout << "error: couldn't write warped peaks into file "
                          << outfile_name << std::endl;
                return -1;
            }
        } else {
            std::cout << "error: unknown file format for file " << input_file
                      << std::endl;
            print_usage();
            return -1;
        }

        // std::vector<double> data;
        // if ((options.find("-parallel") != options.end()) &&
        //(options["-parallel"] == "true" ||
        // options["-parallel"].empty())) {
        // data = Grid::Runners::Parallel::run(max_threads, parameters,
        // all_points);
        //} else {
        // data = Grid::Runners::Serial::run(parameters, all_points);
        //}

        // std::cout << "Saving grid into dat file..." << std::endl;
        // if (!Grid::Files::Dat::write(datfile_stream, data, parameters)) {
        // std::cout << "error: the grid could not be saved properly"
        //<< std::endl;
        // return -1;
        //}
    }

    return 0;
}
