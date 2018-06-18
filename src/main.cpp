#include <iostream>
#include <map>

void print_usage() {
    std::cout << "usage: grid [options] <files>" << std::endl;
}

// TODO: pair<string, string>:  <Flag name, description>. OR:
//  std::tuple<string, string, bool> <Flag name, description, takes_parameters>
const std::map<std::string, std::pair<std::string, bool>> accepted_flags = {
    // Dimensions.
    {"--num_mz", {"*description--", true}},
    {"--num_rt", {"*description--", true}},
    // Bounds.
    {"--min_rt", {"*description--", true}},
    {"--max_rt", {"*description--", true}},
    {"--min_mz", {"*description--", true}},
    {"--max_mz", {"*description--", true}},
    // SmoothingParams.
    {"--smooth_mz", {"*description--", true}},
    {"--sigma_mz", {"*description--", true}},
    {"--sigma_rt", {"*description--", true}},
    // Instrument::Type.
    {"--instrument", {"*description--", true}},
    // Flags.
    {"--warped", {"*description--", false}},
    // Command parameters.
    {"--out_dir", {"*description--", true}},
    {"--help", {"*description--", false}},
};

int main(int argc, char* argv[]) {
    if (argc == 1) {
        print_usage();
        return -1;
    }

    // TODO: Fill argument map.
    // std::vector<std::string> unknown_flags = {};
    for (int i = 1; i < argc; ++i) {
        auto arg = argv[i];
        if (accepted_flags.find(arg) == accepted_flags.end()) {
            std::cout << "unknown option: " << arg << std::endl;
            print_usage();
            return -1;
        }
        auto flag = accepted_flags.at(arg);
        std::string option_arg = "";
        // If the flag takes arguments make sure it's available and fetch it.
        if (flag.second) {
            if (i + 1 >= argc) {
                std::cout << "no parameters specified for " << arg << std::endl;
                print_usage();
                return -1;
            }
            ++i;
            option_arg = argv[i];
            std::cout << arg << " " << option_arg << std::endl;
        }
    }

    // TODO: Check for unknown file format.
    // TODO: Check for no files specified.
    // std::cout << "ARGC: " << argc << std::endl;
    // std::cout << "ARGV:" << std::endl;
    // for (int i = 0; i < argc; ++i) {
    // auto arg = argv[i];
    // std::cout << arg << std::endl;
    //}
    // print_usage();
    return 0;
}
