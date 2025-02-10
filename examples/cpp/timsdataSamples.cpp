#include <iostream>
#include <string>
#include <iomanip>
#include <vector>

// #include "CppSQLite3.h"
// #include "timsdata_cpp.h"
#include "timsdatacpp.h"


using namespace timsdata;


#define OUTPUT_PRECISION 8
#define OUTPUT_WIDTH 9

void dumpFrame(timsdata::TimsData& data, int64_t frameId, int scanEnd);

// When you start this test programs, timsdata.dll is usually not found.
// Either copy the DLL into a place where it is found.
// Or add the following line in the project settings under Debugging/Environment
// PATH=%PATH%;$(SolutionDir)\..\..\win64

// Simple sample program for using the C++ wrapper of the timsdata DLL.
// Please note: this sample is by design very simple, some error handling might be missing.
// Also note: string to/from timsdata are UTF8, this sample does not convert the "native" encoded 
// string to/from UTF-8, it will only work as long as all string only use 7 bit ASCII characters.
int main(int argc, char* argv[])
{
    if(argc != 2)
    {
        std::cerr << "Usage:" << std::endl;
        std::cerr << argv[0] << " <TDF .d directory>" << std::endl;
        return -1;
    }

    std::cout.precision(OUTPUT_PRECISION);

    std::string tdfDirectory(argv[1]);
    // std::string tdfFile = tdfDirectory + "\\analysis.tdf";

    try
    {
        // Open the binary file
        timsdata::TimsData td(tdfDirectory);  // NOTE: UTF-8 conversion needed here!

        // Open the database
        // CppSQLite3DB db;
        // db.open(tdfFile.c_str());  // NOTE: UTF-8 conversion needed here!

        // example: get the number of frames (spectra)
        //auto frames = db.execScalar("SELECT COUNT(*) FROM Frames;");
        auto nframes = td.getNumberOfFrames();
        int frame_id;
        std::string user_input;

        std::cout << "TDF file " << tdfDirectory << " contains " << nframes << " frames." << std::endl;


  
        while (true) {
            std::cout << "Enter a frame ID (must be less than " << nframes << ") or press 'q' to quit: ";
            std::cin >> user_input;

            if (user_input == "q") break;
            
            try {
                frame_id = std::stoi(user_input);
                if (frame_id < 0 || frame_id >= nframes) {
                    std::cerr << "Invalid frame ID. Please enter a value between 0 and " << nframes - 1 << "." << std::endl;
                    continue;
                }
            } catch (const std::exception&) {
                std::cerr << "Invalid input. Please enter a number." << std::endl;
                continue;
            }

            std::string sql_query = "SELECT * FROM Frames WHERE Id = " + std::to_string(frame_id) + ";";
            auto results = td.query(sql_query);

            std::cout << "\nFrame Data:\n";
            for (const auto& row : results) {
                for (const auto& [column, value] : row) {
                    std::cout << std::setw(20) << column << ": ";
                    if (std::holds_alternative<int64_t>(value))
                        std::cout << std::get<int64_t>(value);
                    else if (std::holds_alternative<double>(value))
                        std::cout << std::get<double>(value);
                    else if (std::holds_alternative<std::string>(value))
                        std::cout << std::get<std::string>(value);
                    std::cout << std::endl;
                }
                std::cout << std::string(80, '-') << std::endl;
            }
            dumpFrame(td, frame_id, 200);
        }

        // get information about each frame
        // auto query = db.execQuery("SELECT Id, Time, NumScans FROM Frames;");

        // while(!query.eof())
        // {
        //     auto time     = query.getFloatField("Time");
        //     auto id       = query.getIntField("Id");
        //     auto numScans = query.getIntField("NumScans");

        //     std::cout << "Frame " << id << " has retention time " << time << " seconds" << std::endl;
        //     # dumpFrame(data, id, numScans);

        //     query.nextRow();
        // }
    }
    catch(const std::exception& e)
    {
        std::cerr << "Exception: " << e.what() << std::endl;
    }
    // catch(const CppSQLite3Exception& e)
    // {
    //     std::cerr << "CppSQLite3Exception: " << e.errorMessage() << std::endl;
    // }

	return 0;
}


void dumpFrame(timsdata::TimsData& data, int64_t frameId, int scanEnd)
{
    auto scans = data.readScans(frameId, 0, scanEnd);

    std::cout << std::left << std::setw(10) << "Scan"
              << std::setw(15) << "1/K0"
              << std::setw(15) << "x (index)"
              << std::setw(15) << "x (m/z)"
              << std::setw(15) << "y (intensity)" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    // dump part of the data - for this sample dump a maximum of 15 scans that contain at least one peak
    // and for each of these scans a maximum of 10 peaks
    auto maxScansWithPeaks = std::min(size_t{ 15 }, scans.getNbrScans());
    auto scansWithPeaks = 0u;
    for(auto scan = 0u;  scansWithPeaks < maxScansWithPeaks && scan < scans.getNbrScans(); ++scan)
    {
        auto maxPeaks = std::min(size_t{ 10 }, scans.getNbrPeaks(scan));
        if (maxPeaks == 0) {
            continue;
        } else {
            scansWithPeaks++;
        }
        std::vector<double> mobility;
        data.scanNumToOneOverK0(frameId, { static_cast<double>(scan) }, mobility);
        // std::cout << "Scan " << scan << " --- 1/K0: " << mobility[0] << std::endl;

        auto xAxis = scans.getScanX(scan);
        auto yAxis = scans.getScanY(scan);

        // std::cout << "x (index)    : ";
        // for(size_t counter = 0; counter < maxPeaks; ++counter)
        // {
        //     std::cout << std::setw(OUTPUT_WIDTH) << xAxis.first[counter] << " ";
        // }
        // std::cout << std::endl;

        // std::cout << "x (m/z)      : ";
        std::vector<double> xAxisMasses;
        std::vector<double> indices(xAxis.first, xAxis.second);

        data.indexToMz(frameId, indices, xAxisMasses);

        for (size_t counter = 0; counter < maxPeaks; ++counter) {
            if (counter == 0) {
                std::cout << std::left << std::setw(10) << scan
                          << std::setw(15) << mobility[0];
            } else {
                std::cout << std::setw(25) << ""; // Align columns
            }
            std::cout << std::setw(15) << xAxis.first[counter]
                      << std::setw(15) << xAxisMasses[counter]
                      << std::setw(15) << yAxis.first[counter] << std::endl;
        }

        // for(size_t counter = 0; counter < maxPeaks; ++counter)
        // {
        //     std::cout << std::setw(OUTPUT_WIDTH) << xAxisMasses[counter] << " ";
        // }
        // std::cout << std::endl;

        // // std::cout << "y (intensity): ";
        // for(size_t counter = 0; counter < maxPeaks; ++counter)
        // {
        //     std::cout << std::setw(OUTPUT_WIDTH) << yAxis.first[counter] << " ";
        // }
        // std::cout << std::endl << std::endl;

    }
    std::cout << std::string(70, '-') << std::endl;
}
