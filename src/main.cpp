#include <iostream>
#include <map>
#include <vector>

#include "grid.hpp"

void print_usage() {
    std::cout << "USAGE: grid [-help] [-option x] <files>" << std::endl;
}

// Flag format is map where the key is the flag name and contains a tuple with
// the description and if it takes extra parameters or not: <description,
// takes_parameters>
const std::map<std::string, std::pair<std::string, bool>> accepted_flags = {
    // Dimensions.
    {"-num_mz", {"The number of sampling points for the grid on mz", true}},
    {"-num_rt", {"The number of sampling points for the grid on rt", true}},
    // Bounds.
    {"-min_rt", {"The minimum rt value", true}},
    {"-max_rt", {"The maximum rt value", true}},
    {"-min_mz", {"The minimum mz value", true}},
    {"-max_mz", {"The maximum mz value", true}},
    // SmoothingParams.
    {"-smooth_mz", {"The mass at which the smoothing sigma is given", true}},
    {"-sigma_mz", {"The smoothing sigma in the mz direction", true}},
    {"-sigma_rt", {"The smoothing sigma in the rt direction", true}},
    // Instrument::Type.
    {"-instrument", {"The instrument in which the data was extracted", true}},
    // Flags.
    {"-warped", {"Specify if the output grid will be warped", false}},
    // Command parameters.
    {"-out_dir", {"The output directory", true}},
    {"-help", {"Display available options", false}},
    {"-config", {"Specify the configuration file", true}},
};

bool is_unsigned_integer(std::string& s) { return true; }
bool is_integer(std::string& s) { return true; }
bool is_double(std::string& s) { return true; }

int main(int argc, char* argv[]) {
    if (argc == 1) {
        print_usage();
        return -1;
    }

    // Parse arguments and extract options and file names.
    std::map<std::string, std::string> options;
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
            files.push_back(arg);
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
            for (int i = 0; i < (padding - e.first.size()) + 4; ++i) {
                std::cout << " ";
            }
            std::cout << e.second.first << std::endl;
        }
        return 0;
    }

    if (files.size() == 0) {
        std::cout << "No input files specified." << std::endl;
        print_usage();
        return -1;
    }

    // TODO(alex): If config file is provided, read it and parse it. The
    // parameters specified as command line arguments will override the config
    // file.

    // Parse the options to build the Grid::Parameters struct.
    Grid::Parameters parameters;
    // Get the dimensions.
    // TODO: Accept also -delta_mz/-delta_rt instead of -num_mz/-num_rt
    if ((options.find("-num_mz") == options.end()) ||
        (options.find("-num_rt") == options.end())) {
        std::cout << "Grid dimensions (num_mz, num_rt) not specified"
                  << std::endl;
        return -1;
    }
    auto num_mz = options["-num_mz"];
    auto num_rt = options["-num_rt"];
    if (!is_unsigned_integer(num_mz)) {
        std::cout << "error: "
                  << "num_mz"
                  << " has to be a positive integer" << std::endl;
        print_usage();
        return -1;
    }
    if (!is_unsigned_integer(num_rt)) {
        std::cout << "error: "
                  << "num_rt"
                  << " has to be a positive integer" << std::endl;
        print_usage();
        return -1;
    }
    parameters.dimensions.n = std::stoi(num_mz);
    parameters.dimensions.m = std::stoi(num_rt);

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
    if (!is_double(min_rt)) {
        std::cout << "error: "
                  << "min_rt"
                  << " has to be a number" << std::endl;
        print_usage();
        return -1;
    }
    if (!is_double(max_rt)) {
        std::cout << "error: "
                  << "max_rt"
                  << " has to be a number" << std::endl;
        print_usage();
        return -1;
    }
    if (!is_double(min_mz)) {
        std::cout << "error: "
                  << "min_mz"
                  << " has to be a number" << std::endl;
        print_usage();
        return -1;
    }
    if (!is_double(max_mz)) {
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
    if (!is_double(smooth_mz)) {
        std::cout << "error: "
                  << "smooth_mz"
                  << " has to be a number" << std::endl;
        print_usage();
        return -1;
    }
    if (!is_double(sigma_mz)) {
        std::cout << "error: "
                  << "sigma_mz"
                  << " has to be a number" << std::endl;
        print_usage();
        return -1;
    }
    if (!is_double(sigma_rt)) {
        std::cout << "error: "
                  << "sigma_rt"
                  << " has to be a number" << std::endl;
        print_usage();
        return -1;
    }
    parameters.smoothing_params.mz = std::stod(smooth_mz);
    parameters.smoothing_params.sigma_mz = std::stod(sigma_mz);
    parameters.smoothing_params.sigma_rt = std::stod(min_mz);

    // Get the instrument type.
    if (options.find("-instrument") == options.end()) {
        std::cout << "Instrument type (instrument) not specified" << std::endl;
        print_usage();
        return -1;
    }
    auto instrument = options["-instrument"];
    // Transform instrument to lowercase to prevent typos.
    for (int i = 0; i < instrument.size(); ++i) {
        instrument[i] = std::tolower(instrument[i]);
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

    std::cout << "PRINTING ARGUMENTS:" << std::endl;
    for (const auto& e : options) {
        std::cout << e.first << " " << e.second << std::endl;
    }
    std::cout << "PRINTING FILES:" << std::endl;
    for (const auto& e : files) {
        std::cout << e << std::endl;
    }

    // TODO(alex): Check for unknown file format.
    // TODO(alex): Check for no files specified.
    // TODO(alex): check for file not found.
    return 0;
}
