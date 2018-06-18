#include <iostream>
#include <map>
#include <vector>

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

    // Print help.
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
