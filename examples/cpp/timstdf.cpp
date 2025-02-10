#include "timsdatacpp.h"
#include <iostream>
#include <iomanip>  // for pretty printing
#include <variant>
#include <map>
#include <vector>

using namespace timsdata;

// Helper function to pretty print a map (for the query results)
void pretty_print_row(const std::map<std::string, std::variant<int64_t, double, std::string>>& row) {
    for (const auto& [key, value] : row) {
        std::cout << std::setw(20) << std::left << key << ": ";
        
        // Print the variant depending on its type
        std::visit([](auto&& arg) { std::cout << arg << std::endl; }, value);
    }
}

// Main function
int main(int argc, char* argv[]) {
    // Check if the required argument (analysis directory) is passed
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <path_to_.d_folder> [sql_query]\n";
        return 1;
    }

    // Get the .d folder path from the command-line argument
    std::string analysis_directory_name = argv[1];

    // Check if an SQL query is provided
    std::string sql_query;
    if (argc >= 3) {
        sql_query = argv[2];  // Use the provided query
    }

    // Initialize options for TimsData
    bool use_recalibration = false; // Default value
    pressure_compensation_strategy pressure_compensation = AnalyisGlobalPressureCompensation;  // Default

    try {
        // Initialize the TimsData object with the provided .d folder
        TimsData tims_data(analysis_directory_name, use_recalibration, pressure_compensation);

        // 1. Get and print the list of tables
        std::cout << "Tables in the database:\n";
        std::vector<std::string> tables = tims_data.get_tables();
        for (const auto& table : tables) {
            std::cout << "  " << table << std::endl;
        }

        // 2. Get and print the schema of each table
        std::cout << "\nTable schemas:\n";
        for (const auto& table : tables) {
            std::cout << "\nSchema for table: " << table << std::endl;
            std::string schema = tims_data.get_schema(table);
            std::cout << schema << std::endl;
        }

        // 3. If an SQL query was provided, execute and print results
        if (!sql_query.empty()) {
            std::cout << "\nExecuting query: " << sql_query << "\n";
            std::vector<std::map<std::string, std::variant<int64_t, double, std::string>>> query_results = tims_data.query(sql_query);

            std::cout << "\nQuery Results:\n";
            for (const auto& row : query_results) {
                pretty_print_row(row);
                std::cout << "-----------------------------\n";
            }
        }

        tims_data.relationships();

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

