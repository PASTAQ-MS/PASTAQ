#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <thread>
#include <vector>

#include "grid.hpp"
#include "grid_files.hpp"
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

void parse_json(const std::filesystem::path& path, options_map& options,
                std::vector<std::string>& files) {
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
    std::regex grid_regex("\"grid\"[[:space:]]*:[[:space:]]*\\{(.*)\\}");
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
            auto begin = content.find("[", pos) + 1;
            auto end = content.find("]", begin) - 1;
            if (begin == std::string::npos || end == std::string::npos ||
                begin < 0 || end - begin < 0) {
                std::cout << "error: malformed config file" << std::endl;
                std::exit(-1);
            }

            auto config_files = content.substr(begin, end - begin);
            for (auto& ch : config_files) {
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
            auto begin = content.find("{", pos) + 1;
            auto end = content.find("}", begin) - 1;
            if (begin == std::string::npos || end == std::string::npos ||
                begin < 0 || end - begin < 0) {
                std::cout << "error: malformed config file" << std::endl;
                std::exit(-1);
            }

            auto config_parameters = content.substr(begin, end - begin);
            for (auto& ch : config_parameters) {
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
                name = "-" + name;
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
            auto begin = content.find("{", pos) + 1;
            auto end = content.find("}", begin) - 1;
            if (begin == std::string::npos || end == std::string::npos ||
                begin < 0 || end - begin < 0) {
                std::cout << "error: malformed config file" << std::endl;
                std::exit(-1);
            }

            auto config_config = content.substr(begin, end - begin);
            for (auto& ch : config_config) {
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
                name = "-" + name;
                if (options.find(name) == options.end()) {
                    options[name] = value;
                }
            }
        }
    }
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

// Splits the parameters into n_split sections of the same dimensions.n.
std::vector<Grid::Parameters> split_parameters(
    const Grid::Parameters& original_params, unsigned int n_splits) {
    // In order to determine the overlapping of the splits we need to calculate
    // what is the maximum distance that will be used by the kernel smoothing.
    // To avoid aliasing we will overlap at least the maximum kernel width.
    //
    // The kernel in rt is always 2 * sigma_rt in both directions.
    unsigned int kernel_width = Grid::y_index(
        original_params.bounds.min_rt + 4 * Grid::sigma_rt(original_params),
        original_params);

    // We need to make sure that we have the minimum number of points for the
    // splits. Since we have an overlap of a single kernel_width, we need to
    // have at least twice that amount of points in order to support full
    // overlap in both directions.
    unsigned int min_segment_width = 2 * kernel_width;
    unsigned int segment_width = original_params.dimensions.m / n_splits;
    if (segment_width < min_segment_width) {
        segment_width = min_segment_width;
    }

    // If the orginal parameters don't contain the required minimum number of
    // points for segmentation, we can only use one segment.
    if ((original_params.dimensions.m - 1) < min_segment_width) {
        return std::vector<Grid::Parameters>({original_params});
    }

    // How many segments do we have with the given segment_width.
    unsigned int num_segments = original_params.dimensions.m / segment_width;
    if (original_params.dimensions.m % segment_width) {
        // If we need more segments that the maximum we specify we need to try
        // one less split and adjust the sizes accordingly.
        if (num_segments + 1 > n_splits) {
            segment_width = original_params.dimensions.m / (n_splits - 1);
            num_segments = original_params.dimensions.m / segment_width;
            if (original_params.dimensions.m % segment_width) {
                ++num_segments;
            }
        }
    }

    std::vector<Grid::Parameters> all_parameters;
    auto min_i = 0;
    auto max_i = 0;
    for (size_t i = 0; i < num_segments; ++i) {
        if (i == 0) {
            min_i = 0;
        } else {
            min_i = segment_width * i - kernel_width;
        }
        max_i = segment_width * (i + 1) - 1;
        if (max_i > original_params.dimensions.m) {
            max_i = original_params.dimensions.m - 1;
        }
        // Prepare the next Grid::Parameters object.
        Grid::Parameters parameters(original_params);
        parameters.bounds.min_rt = Grid::rt_at(min_i, original_params).value();
        parameters.bounds.max_rt = Grid::rt_at(max_i, original_params).value();
        parameters.dimensions.m = max_i - min_i + 1;
        all_parameters.emplace_back(parameters);
    }
    return all_parameters;
}

// TODO(alex): this is very memory inefficient, we should avoid copying the
// array of peaks into the different groups. We can either store the references
// or move the values.
std::vector<std::vector<Grid::Peak>> assign_peaks(
    const std::vector<Grid::Parameters>& all_parameters,
    const std::vector<Grid::Peak>& peaks) {
    std::vector<std::vector<Grid::Peak>> groups(all_parameters.size());

    for (const auto& peak : peaks) {
        for (size_t i = 0; i < all_parameters.size(); ++i) {
            auto parameters = all_parameters[i];
            double sigma_rt = Grid::sigma_rt(parameters);
            if (peak.rt + 2 * sigma_rt < parameters.bounds.max_rt) {
                groups[i].push_back(peak);
                break;
            }
            // If we ran out of parameters this peak is assigned to the last
            // one.
            if (i == all_parameters.size() - 1) {
                groups[i].push_back(peak);
            }
        }
    }
    return groups;
}

std::vector<double> merge_groups(
    std::vector<Grid::Parameters>& parameters_array,
    std::vector<std::vector<double>>& data_array) {
    std::vector<double> merged;
    // Early return if there are errors.
    // TODO(alex): we could return a std::nullopt here.
    if (data_array.empty() || parameters_array.empty() ||
        parameters_array.size() != data_array.size()) {
        return merged;
    }
    merged.insert(end(merged), begin(data_array[0]), end(data_array[0]));

    for (size_t n = 1; n < data_array.size(); ++n) {
        auto beg_next = 1 + Grid::y_index(parameters_array[n - 1].bounds.max_rt,
                                          parameters_array[n]);
        // Sum the overlapping sections.
        int i = merged.size() - beg_next * parameters_array[n - 1].dimensions.n;
        for (size_t j = 0;
             j < (parameters_array[n - 1].dimensions.n * beg_next); ++j) {
            merged[i] += data_array[n][j];
            ++i;
        }
        // Insert the next slice.
        merged.insert(
            end(merged),
            begin(data_array[n]) + beg_next * parameters_array[n].dimensions.n,
            end(data_array[n]));
    }
    return merged;
}

std::vector<double> run_parallel(unsigned int max_threads,
                                 const Grid::Parameters& parameters,
                                 const std::vector<Grid::Peak>& all_peaks) {
    // Split parameters and peaks into the corresponding  groups.
    auto all_parameters = split_parameters(parameters, max_threads);
    auto groups = assign_peaks(all_parameters, all_peaks);
    std::cout << "Indexes size: " << groups.size() << std::endl;
    if (groups.size() != all_parameters.size()) {
        std::cout << "error: could not divide the peaks into " << max_threads
                  << " groups" << std::endl;
        std::exit(-1);
    }

    // Allocate data memory.
    std::cout << "Allocating memory..." << std::endl;
    std::vector<std::vector<double>> data_array;
    for (const auto& parameters : all_parameters) {
        data_array.emplace_back(std::vector<double>(parameters.dimensions.n *
                                                    parameters.dimensions.m));
    }

    // Splatting!
    std::cout << "Splatting peaks into concurrent groups..." << std::endl;
    std::vector<std::thread> threads;
    for (size_t i = 0; i < groups.size(); ++i) {
        threads.push_back(
            std::thread([&groups, &all_parameters, &data_array, i]() {
                // Perform splatting in this group.
                for (const auto& peak : groups[i]) {
                    Grid::splat(peak, all_parameters[i], data_array[i]);
                }
            }));
    }

    // Wait for the threads to finish.
    for (auto& thread : threads) {
        thread.join();
    }

    std::cout << "Merging concurrent groups..." << std::endl;
    return merge_groups(all_parameters, data_array);
}

std::vector<double> run_serial(const Grid::Parameters& parameters,
                               const std::vector<Grid::Peak>& all_peaks) {
    // Instantiate memory.
    std::vector<double> data(parameters.dimensions.n * parameters.dimensions.m);

    // Perform grid splatting.
    std::cout << "Splatting peaks into grid..." << std::endl;
    for (const auto& peak : all_peaks) {
        Grid::splat(peak, parameters, data);
    }
    return data;
}

int main(int argc, char* argv[]) {
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

    if (options.find("-warped") != options.end() &&
        (options["-warped"] == "true" || options["-warped"] == "")) {
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
         (options["-parallel"] == "true" || options["-parallel"] == ""))) {
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
    for (const auto& file_name : files) {
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
        for (auto& ch : lowercase_extension) {
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

            // Prepare the name of the rawdump output file.
            // TODO(alex): this should be optional.
            auto rawdump_name =
                options["-out_dir"] /
                input_file.filename().replace_extension(".rawdump");
            std::ofstream rawdump_stream;
            rawdump_stream.open(rawdump_name, std::ios::out | std::ios::binary);
            if (!rawdump_stream) {
                std::cout << "error: could not open file " << rawdump_name
                          << " for writing" << std::endl;
                return -1;
            }

            std::cout << "Parsing file..." << std::endl;
            auto peaks = XmlReader::read_next_scan(stream, parameters);
            if (peaks == std::nullopt) {
                std::cout << "error: no peaks found on file " << input_file
                          << " for the given parameters" << std::endl;
                return -1;
            }
            std::vector<Grid::Peak> all_peaks = {};
            do {
                all_peaks.insert(end(all_peaks), begin(peaks.value()),
                                 end(peaks.value()));
                peaks = XmlReader::read_next_scan(stream, parameters);
            } while (peaks != std::nullopt);

            // TODO(alex): this should be optional.
            if (!Grid::Files::Rawdump::write(rawdump_stream, all_peaks)) {
                std::cout << "error: the raw dump could not be saved properly"
                          << std::endl;
                return -1;
            }

            std::vector<double> data;
            if ((options.find("-parallel") != options.end()) &&
                (options["-parallel"] == "true" ||
                 options["-parallel"] == "")) {
                data = run_parallel(max_threads, parameters, all_peaks);
            } else {
                data = run_serial(parameters, all_peaks);
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

            // Load the peaks into memory.
            std::vector<Grid::Peak> all_peaks;
            if (!Grid::Files::Rawdump::read(stream, all_peaks)) {
                std::cout << "error: the raw dump could not be loaded properly"
                          << std::endl;
                return -1;
            }
            if (all_peaks.empty()) {
                std::cout << "error: the raw dump does not contain any peaks"
                          << std::endl;
                return -1;
            }
            std::cout << "Loaded " << all_peaks.size() << " peaks" << std::endl;

            std::vector<double> data;
            if ((options.find("-parallel") != options.end()) &&
                (options["-parallel"] == "true" ||
                 options["-parallel"] == "")) {
                data = run_parallel(max_threads, parameters, all_peaks);
            } else {
                data = run_serial(parameters, all_peaks);
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
