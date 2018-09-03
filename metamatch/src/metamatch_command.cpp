#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <thread>
#include <vector>

#include "metamatch/metamatch.hpp"

// Type aliases.
using options_map = std::map<std::string, std::string>;

void print_usage() {
    std::cout << "USAGE: metamatch [-help] [options]" << std::endl;
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

    // Find the contents of the configuration.
    std::regex config_regex(R"("metamatch"[[:space:]]*:[[:space:]]*\{(.*)\})");
    std::smatch matches;
    std::regex_search(content, matches, config_regex);
    if (matches.size() != 2 || matches[1].str().empty()) {
        std::cout << "error: could not find \"metamatch\" on the config file"
                  << std::endl;
        std::exit(-1);
    }
    content = matches[1];

    // Match the different parts of the configuration. Note that we are
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
        // MetaMatch parameters.
        {"-file_list",
         {"A file containing the paths and classes of the aligned peak lists",
          true}},
        {"-radius_mz",
         {"The maximum distance in mz that can be used for clustering", true}},
        {"-radius_rt",
         {"The maximum distance in rt that can be used for clustering", true}},
        {"-fraction",
         {"The percentage of samples that must contain non zero values to "
          "consider a cluster valid [0.0,1.0]",
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

    // Parse the options to build the Grid::Parameters struct.
    MetaMatch::Parameters parameters = {};

    // Get the MetaMatch parameters.
    if (options.find("-radius_mz") == options.end()) {
        std::cout << "Radius mz (radius_mz) not specified" << std::endl;
        return -1;
    }
    auto radius_mz = options["-radius_mz"];
    if (!is_number(radius_mz)) {
        std::cout << "error: radius_mz has to be a number" << std::endl;
        print_usage();
        return -1;
    }

    parameters.radius_mz = std::stod(radius_mz);
    if (options.find("-radius_rt") == options.end()) {
        std::cout << "Radius rt (radius_rt) not specified" << std::endl;
        return -1;
    }
    auto radius_rt = options["-radius_rt"];
    if (!is_number(radius_rt)) {
        std::cout << "error: radius_rt has to be a number" << std::endl;
        print_usage();
        return -1;
    }
    parameters.radius_rt = std::stod(radius_rt);

    if (options.find("-fraction") == options.end()) {
        std::cout << "Radius rt (fraction) not specified" << std::endl;
        return -1;
    }
    auto fraction = options["-fraction"];
    if (!is_number(fraction)) {
        std::cout << "error: fraction has to be a number" << std::endl;
        print_usage();
        return -1;
    }
    parameters.fraction = std::stod(fraction);
    if (parameters.fraction < 0 || parameters.fraction > 1) {
        std::cout << "error: fraction has to be a number between 0 and 1"
                  << std::endl;
        print_usage();
        return -1;
    }

    if (options.find("-file_list") == options.end()) {
        std::cout << "File list (file_list) not specified" << std::endl;
        return -1;
    }
    auto file_list = options["-file_list"];

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

    // Read the file list.
    std::cout << "Reading file list: " << file_list << std::endl;
    {
        std::ifstream stream;
        stream.open(file_list);
        if (!stream) {
            std::cout << "error: could not open file list " << file_list
                      << std::endl;
            return -1;
        }

        // Prepare the name of the output file.
        // TODO(alex): Only check if directory exists, do not open the file for
        // writing until we are done with the process.
        std::filesystem::path output_file_name = "metapeaks.csv";
        auto outfile_name = options["-out_dir"] / output_file_name;
        std::ofstream outfile_stream;
        outfile_stream.open(outfile_name, std::ios::out | std::ios::binary);
        if (!outfile_stream) {
            std::cout << "error: could not open file " << outfile_name
                      << " for writing" << std::endl;
            return -1;
        }

        auto files = MetaMatch::read_file_list(stream);
        std::vector<MetaMatch::Peak> metapeaks;
        size_t file_id = 0;
        std::vector<size_t> classes;
        for (const auto &[file, class_id] : files) {
            std::filesystem::path input_file = file;
            std::cout << "Reading file: " << input_file << std::endl;
            std::ifstream peaks_stream;
            peaks_stream.open(input_file);
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
                std::cout << "Reading peaks from file: " << input_file
                          << std::endl;
                if (!Centroid::Files::Bpks::read_peaks(peaks_stream,
                                                       &grid_params, &peaks)) {
                    std::cout << "error: couldn't read peaks from the file list"
                              << std::endl;
                    return -1;
                }
            } else if (lowercase_extension == ".csv") {
                // TODO: read into peak array.
                if (!Centroid::Files::Csv::read_peaks(peaks_stream, &peaks)) {
                    std::cout << "error: couldn't read peaks from the file list"
                              << std::endl;
                    return -1;
                }
            } else {
                std::cout << "error: unknown file format for file "
                          << input_file << std::endl;
                print_usage();
                return -1;
            }
            for (const auto &peak : peaks) {
                metapeaks.push_back(
                    {peak, file_id, class_id, -1, peak.mz, peak.rt});
            }
            ++file_id;

            if (std::find(classes.begin(), classes.end(), class_id) ==
                classes.end()) {
                classes.push_back(class_id);
            }
        }

        // Execute MetaMatch here.
        // TODO(alex): Error checking!
        std::cout << "Finding candidates..." << std::endl;
        MetaMatch::find_candidates(metapeaks, parameters);
        // TODO(alex): Error checking!
        auto orphans = MetaMatch::extract_orphans(metapeaks);
        std::cout << "Extracting orphans..." << std::endl;
        // TODO(alex): Error checking!
        std::cout << "Building cluster table..." << std::endl;
        auto clusters = MetaMatch::reduce_cluster(metapeaks, files.size());
        // TODO(alex): Error checking!
        std::cout << "Writing table to disk..." << std::endl;
        MetaMatch::write_clusters(outfile_stream, clusters, files.size());
    }

    return 0;
}
