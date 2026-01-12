#include "timsdatacpp.h"
#include <iostream>
#include <iomanip>  // for pretty printing
#include <variant>
#include <map>
#include <vector>

using namespace timsdata;

// Helper function to pretty print a frame (assuming Frames have ID and some metadata)
// void pretty_print_frames(const std::vector<std::map<std::string, std::variant<int64_t, double, std::string>>>& frames) {
//     for (const auto& frame : frames) {
//         for (const auto& [key, value] : frame) {
//             std::cout << std::setw(20) << std::left << key << ": ";
//             std::visit([](auto&& arg) { std::cout << arg << std::endl; }, value);
//         }
//         std::cout << "-----------------------------\n";
//     }
// }

// Main function
int main(int argc, char* argv[]) {
    // Check if the required argument (analysis directory) is passed
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <path_to_.d_folder>\n";
        return 1;
    }

    // Get the .d folder path from the command-line argument
    std::string analysis_directory_name = argv[1];

    // Initialize options for TimsData
    bool use_recalibration = false; // Default value
    pressure_compensation_strategy pressure_compensation = AnalyisGlobalPressureCompensation;  // Default

    try {
        // Initialize the TimsData object with the provided .d folder
        TimsData tims_data(analysis_directory_name, use_recalibration, pressure_compensation);

        // Print tables for reference
        std::cout << "Tables in the database:\n";
        std::vector<std::string> tables = tims_data.get_tables();
        for (const auto& table : tables) {
            std::cout << "  " << table << std::endl;
        }

        // Test populateFrames with the given filter
        std::cout << "\nTesting populateFrames(\"MsMsType = 0\")\n";
        tims_data.populateFrames("MsMsType = 0");

        // Print the number of frames after filtering
        std::cout << "Number of frames loaded: " << tims_data.getNumberOfFrames() << std::endl;

        // // Print first few frames (if possible) for verification
        // std::vector<std::map<std::string, std::variant<int64_t, double, std::string>>> frames = tims_data.getFramesData();
        // if (!frames.empty()) {
        //     std::cout << "\nFirst few frames:\n";
        //     pretty_print_frames(frames);
        // } else {
        //     std::cout << "No frames matched the filter.\n";
        // }

        // tims_data.relationships();

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}