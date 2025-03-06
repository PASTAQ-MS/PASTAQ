#include <stdexcept>
#include <string>
#include <cstdint>
#include <numeric>
#include <vector>
#include <limits>
#include <memory>
#include <utility> // For std::move
#include <iostream>
#include <iomanip>

// #include <cstddef>
// #include "CppSQLite3.h"

// #include "timsdata.h"
#include "timsdatacpp.h"


namespace timsdata
{
    size_t FrameProxy::getNbrScans() const {
        return num_scans;
    }

    size_t FrameProxy::getTotalNbrPeaks() const {
        return scan_offsets.back();
    }

    size_t FrameProxy::getNbrPeaks(size_t scan_num) const {
        throwIfInvalidScanNumber(scan_num);
        return pData[scan_num];
    }

    FrameProxy::FrameIteratorRange FrameProxy::getScanX(size_t scan_num) const {
        return makeRange(scan_num, 0);
    }

    FrameProxy::FrameIteratorRange FrameProxy::getScanY(size_t scan_num) const {
        return makeRange(scan_num, pData[scan_num]);
    }

    void FrameProxy::throwIfInvalidScanNumber(size_t scan_num) const {
        if(scan_num >= getNbrScans())
            throw std::invalid_argument("Scan number out of range.");
    }

    FrameProxy::FrameIteratorRange FrameProxy::makeRange(size_t scan_num, size_t offset) const {
        throwIfInvalidScanNumber(scan_num);
        const uint32_t* p = pData.get() + num_scans + 2 * scan_offsets[scan_num] + offset;
        return std::make_pair(p, p + pData[scan_num]);
    }

    std::string getLastError()
    {
        uint32_t len = tims_get_last_error_string(0, 0);
        std::unique_ptr<char[]> buf(new char[len]);
        tims_get_last_error_string(buf.get(), len);
        return std::string(buf.get(), buf.get() + len - 1);
    }

    void throwLastError()
    {
        throw std::runtime_error(getLastError());
    }

    TimsData::TimsData(const std::string& analysis_directory_name, bool use_recalibration, pressure_compensation_strategy pressure_compensation)
        : handle(0)
        , initial_frame_buffer_size(128)
        , tdfFile(analysis_directory_name + "/analysis.tdf")
    {
        handle = tims_open_v2(analysis_directory_name.c_str(), use_recalibration, pressure_compensation);
        if(handle == 0)
            throwLastError();

        // Initialize the SQLite database connection
        try {
            db.open(tdfFile.c_str()); // Opens the database in read mode
        } catch (const CppSQLite3Exception& e) {
            // Handle the error if the database connection fails
            std::cerr << "Failed to open SQLite database: " << e.what() << std::endl;
            throw; // Re-throw the exception after logging it
        }
    }

    TimsData::~TimsData()
    {
        tims_close(handle);
        db.close();
    }

    std::vector<std::string> TimsData::get_tables() const {
        std::vector<std::string> tables;
        CppSQLite3Query query = db.execQuery("SELECT name FROM sqlite_master WHERE type='table';");
        while (!query.eof()) {
            tables.push_back(query.getStringField(0));
            query.nextRow();
        }
        return tables;
    }

    std::string TimsData::get_schema(const std::string& table_name) const {
        std::ostringstream schema;
        std::string query_str = "PRAGMA table_info(" + table_name + ");";  // Create a string
        CppSQLite3Query query = db.execQuery(query_str.c_str());

        schema << "Column name" << std::setw(21) << "Data type\n";
        schema << std::string(30, '-') << std::endl;
        while (!query.eof()) {
            size_t length = std::string(query.getStringField(1)).length();
            schema << query.getStringField(1) << std::setw(23-length) << " " << query.getStringField(2) << std::endl;
            query.nextRow();
        }
        schema << std::string(30, '-') << std::endl;
        return schema.str();
    }

    std::vector<std::map<std::string, std::variant<int64_t, double, std::string>>> TimsData::query(const std::string& sql_query) const {
        std::vector<std::map<std::string, std::variant<int64_t, double, std::string>>> results;
        CppSQLite3Query query = db.execQuery(sql_query.c_str());
        int col_count = query.numFields();
        while (!query.eof()) {
            std::map<std::string, std::variant<int64_t, double, std::string>> row;
            for (int i = 0; i < col_count; ++i) {
                std::string col_name = query.fieldName(i);
                if (query.fieldDataType(i) == SQLITE_INTEGER) {
                    row[col_name] = static_cast<int64_t>(query.getInt64Field(i));
                } else if (query.fieldDataType(i) == SQLITE_FLOAT) {
                    row[col_name] = query.getFloatField(i);
                } else {
                    row[col_name] = query.getStringField(i);
                }
            }
            results.push_back(row);
            query.nextRow();
        }
        return results;
    }

    uint64_t TimsData::getHandle() const
    {
        return handle;
    }

    FrameProxy TimsData::readScans(int64_t frame_id, uint32_t scan_begin, uint32_t scan_end)
    {
        if(scan_end < scan_begin)
            throw std::runtime_error("scan_end must be >= scan_begin");

        const uint32_t num_scans = scan_end - scan_begin;
        std::unique_ptr<uint32_t[]> pData;

        for(;;) {
            pData.reset(new uint32_t[initial_frame_buffer_size]);

            uint32_t required_len = tims_read_scans_v2(handle, frame_id, scan_begin, scan_end, pData.get(), uint32_t(4 * initial_frame_buffer_size));
            if(required_len == 0)
                throwLastError();

            if(4 * initial_frame_buffer_size > required_len) {
                if(required_len < 4 * num_scans)
                    throw std::runtime_error("Data array too small.");
                return FrameProxy(num_scans, std::move(pData));
            }

            if(required_len > 16777216)
                throw std::runtime_error("Maximum expected frame size exceeded.");

            initial_frame_buffer_size = required_len / 4 + 1;
        }
    }

    void TimsData::doTransformation(int64_t frame_id, const std::vector<double>& in, std::vector<double>& out, BdalTimsConversionFunction* func)
    {
        if(in.empty())
        {
            out.clear();
            return;
        }
        if(in.size() > std::numeric_limits<uint32_t>::max())
            throw std::runtime_error("Input range too large.");
        out.resize(in.size());
        func(handle, frame_id, &in[0], &out[0], uint32_t(in.size()));
    }

    /// Returns the number of frames (spectra) in the TIMS analysis.
    ///
    /// \throws std::exception if the database query fails
    uint32_t TimsData::getNumberOfFrames() const
    {
        // Query the number of frames
        auto number_of_frames = db.execScalar("SELECT COUNT(*) FROM Frames;");
        return static_cast<uint32_t>(number_of_frames);
    }


    std::pair<std::vector<std::string>, std::vector<std::vector<std::variant<int64_t, double, std::string>>>> 
TimsData::getFramesTable() const {
    // std::vector<std::vector<std::variant<int64_t, double, std::string>>> TimsData::getFramesTable() const {
        std::vector<std::string> column_names;
        std::vector<std::vector<std::variant<int64_t, double, std::string>>> frames_data;

        try {
            // Execute SQL query to fetch all rows from the Frames table
            CppSQLite3Query query = db.execQuery("SELECT * FROM Frames;");

            // Retrieve column names
            int num_fields = query.numFields();
            for (int col = 0; col < num_fields; ++col) {
                column_names.push_back(query.fieldName(col));
            }

            // Loop through each row in the query result
            while (!query.eof()) {
                std::vector<std::variant<int64_t, double, std::string>> frame_row;

                // Loop through each column in the current row
                for (int col = 0; col < query.numFields(); ++col) {
                    const char* value = query.fieldValue(col);
                    const char* colType = query.fieldDeclType(col);
    
                    if (!value) {
                        frame_row.push_back("");  // Handle NULL values as empty strings
                    } else if (std::string(colType).find("INT") != std::string::npos) {
                        frame_row.push_back(static_cast<int64_t>(std::stoll(value)));
                    } else if (std::string(colType).find("REAL") != std::string::npos || 
                           std::string(colType).find("FLOAT") != std::string::npos) {
                        frame_row.push_back(std::stod(value));
                    } else {
                        frame_row.push_back(std::string(value));  // Default to string
                    }
                }

                frames_data.push_back(frame_row);
                query.nextRow();
            }
        } catch (const CppSQLite3Exception& e) {
            throw std::runtime_error("Error fetching data from Frames table: " + std::string(e.errorMessage()));
        }

        return {column_names, frames_data};
        // return frames_data;
    }

    void TimsData::relationships() const {
        try {
            // Retrieve the list of tables
            CppSQLite3Query tableQuery = db.execQuery("SELECT name FROM sqlite_master WHERE type='table';");
            std::map<std::string, std::vector<std::pair<std::string, std::string>>> relationships;

            while (!tableQuery.eof()) {
                std::string table_name = tableQuery.getStringField(0);
                std::string query = "PRAGMA foreign_key_list('" + table_name + "');";
                CppSQLite3Query fkQuery = db.execQuery(query.c_str());

                while (!fkQuery.eof()) {
                    std::string source_table = fkQuery.getStringField(2); // Referenced table
                    std::string target_table = table_name; // Table containing the foreign key
                    std::string column_relationship = std::string(fkQuery.getStringField(3)) + " -> " + fkQuery.getStringField(4);
                    relationships[source_table].emplace_back(target_table, column_relationship);
                    fkQuery.nextRow();
                }
                tableQuery.nextRow();
            }

            // Print relationships
            std::cout << "Table Relationships (Foreign Keys):\n" << std::endl;
            for (const auto& entry : relationships) {
                std::cout << "Table '" << entry.first << "' is referenced by:" << std::endl;
                for (const auto& rel : entry.second) {
                    std::cout << "  -> " << rel.first << " (" << rel.second << ")" << std::endl;
                }
                std::cout << std::endl;
            }
        } catch (const CppSQLite3Exception& e) {
            std::cerr << "SQLite error: " << e.errorMessage() << std::endl;
        }
    } // relationships

    void TimsData::populateFrames(const std::string& filter) {
        // Construct the SQL query with optional filter
        std::string query_str = "SELECT * FROM Frames";
        if (!filter.empty()) {
            query_str += " WHERE " + filter;
        }
        query_str += ";";

        std::vector<std::string> column_names;
        std::vector<std::vector<std::variant<int64_t, double, std::string>>> frames_data;

        try {
            // Execute SQL query to fetch frames
            CppSQLite3Query query = db.execQuery(query_str.c_str());

            // Retrieve column names
            int num_fields = query.numFields();
            for (int col = 0; col < num_fields; ++col) {
                column_names.push_back(query.fieldName(col));
            }

            // Loop through each row in the query result
            while (!query.eof()) {
                std::vector<std::variant<int64_t, double, std::string>> frame_row;

                // Loop through each column in the current row
                for (int col = 0; col < query.numFields(); ++col) {
                    const char* value = query.fieldValue(col);
                    const char* colType = query.fieldDeclType(col);

                    if (!value) {
                        frame_row.push_back("");  // Handle NULL values as empty strings
                    } else if (std::string(colType).find("INT") != std::string::npos) {
                        frame_row.push_back(static_cast<int64_t>(std::stoll(value)));
                    } else if (std::string(colType).find("REAL") != std::string::npos ||
                               std::string(colType).find("FLOAT") != std::string::npos) {
                        frame_row.push_back(std::stod(value));
                    } else {
                        frame_row.push_back(std::string(value));  // Default to string
                    }
                }

                frames_data.push_back(frame_row);
                query.nextRow();
            }
        } catch (const CppSQLite3Exception& e) {
            throw std::runtime_error("Error fetching data from Frames table: " + std::string(e.errorMessage()));
        }

        // Map column names to indices for easier reference
        std::map<std::string, size_t> column_map;
        for (size_t i = 0; i < column_names.size(); ++i) {
            column_map[column_names[i]] = i;
        }

        // Clear existing frames in the cache before populating
        frames_.clear();

        // Loop through each row of the frames data
        for (const auto& row : frames_data) {
            Frame frame;

            // Populate Frame fields, checking for column existence
            if (column_map.count("Id") && std::holds_alternative<int64_t>(row[column_map["Id"]]))
                frame.setId(std::get<int64_t>(row[column_map["Id"]]));

            if (column_map.count("Time") && std::holds_alternative<double>(row[column_map["Time"]]))
                frame.setTime(std::get<double>(row[column_map["Time"]]));

            if (column_map.count("Polarity") && std::holds_alternative<std::string>(row[column_map["Polarity"]]))
                frame.setPolarity(std::get<std::string>(row[column_map["Polarity"]]));

            if (column_map.count("ScanMode") && std::holds_alternative<int64_t>(row[column_map["ScanMode"]]))
                frame.setScanMode(static_cast<int>(std::get<int64_t>(row[column_map["ScanMode"]])));

            if (column_map.count("MsMsType") && std::holds_alternative<int64_t>(row[column_map["MsMsType"]]))
                frame.setMsMsType(static_cast<int>(std::get<int64_t>(row[column_map["MsMsType"]])));

            if (column_map.count("TimsId") && std::holds_alternative<int64_t>(row[column_map["TimsId"]]))
                frame.setTimsId(static_cast<int32_t>(std::get<int64_t>(row[column_map["TimsId"]])));

            if (column_map.count("MaxIntensity") && std::holds_alternative<int64_t>(row[column_map["MaxIntensity"]]))
                frame.setMaxIntensity(static_cast<uint32_t>(std::get<int64_t>(row[column_map["MaxIntensity"]])));

            if (column_map.count("SummedIntensities") && std::holds_alternative<int64_t>(row[column_map["SummedIntensities"]]))
                frame.setSummedIntensities(static_cast<uint64_t>(std::get<int64_t>(row[column_map["SummedIntensities"]])));

            if (column_map.count("NumScans") && std::holds_alternative<int64_t>(row[column_map["NumScans"]]))
                frame.setNumScans(static_cast<uint32_t>(std::get<int64_t>(row[column_map["NumScans"]])));

            if (column_map.count("NumPeaks") && std::holds_alternative<int64_t>(row[column_map["NumPeaks"]]))
                frame.setNumPeaks(static_cast<uint32_t>(std::get<int64_t>(row[column_map["NumPeaks"]])));

            if (column_map.count("MzCalibration") && std::holds_alternative<double>(row[column_map["MzCalibration"]]))
                frame.setMzCalibration(std::get<double>(row[column_map["MzCalibration"]]));

            if (column_map.count("T1") && std::holds_alternative<double>(row[column_map["T1"]]))
                frame.setT1(std::get<double>(row[column_map["T1"]]));

            if (column_map.count("T2") && std::holds_alternative<double>(row[column_map["T2"]]))
                frame.setT2(std::get<double>(row[column_map["T2"]]));

            if (column_map.count("TimsCalibration") && std::holds_alternative<double>(row[column_map["TimsCalibration"]]))
                frame.setTimsCalibration(std::get<double>(row[column_map["TimsCalibration"]]));

            if (column_map.count("PropertyGroup") && std::holds_alternative<int64_t>(row[column_map["PropertyGroup"]]))
                frame.setPropertyGroup(static_cast<int>(std::get<int64_t>(row[column_map["PropertyGroup"]])));
            else
                frame.setPropertyGroup(-1); // Use -1 to indicate NULL

            if (column_map.count("AccumulationTime") && std::holds_alternative<double>(row[column_map["AccumulationTime"]]))
                frame.setAccumulationTime(std::get<double>(row[column_map["AccumulationTime"]]));

            if (column_map.count("RampTime") && std::holds_alternative<double>(row[column_map["RampTime"]]))
                frame.setRampTime(std::get<double>(row[column_map["RampTime"]]));

            if (column_map.count("Pressure") && std::holds_alternative<double>(row[column_map["Pressure"]]))
                frame.setPressure(std::get<double>(row[column_map["Pressure"]]));

            // Store the frame in the cache
            frames_[frame.getId()] = std::move(frame);
        }
    }    //populateFrames

    // void TimsData::populateSpectra() {
    //     // Ensure frames are populated first
    //     if (frames_.empty()) {
    //         throw std::runtime_error("Frames are not populated. Call populateFrames() first.");
    //     }

    //     // Loop through all frames
    //     for (auto& frame_pair : frames_) {
    //         Frame& frame = frame_pair.second;

    //         // Clear existing scans to avoid duplication
    //         // frame.clearScans();
    //         frame.clearSpectra();

    //         // Read scans for the current frame
    //         FrameProxy frameProxy = readScans(frame.getId(), 0, frame.getNumScans());

    //         for (uint32_t scan_idx = 0; scan_idx < frameProxy.getNbrScans(); ++scan_idx) {
    //             // Skip scans without peaks
    //             auto numberOfPeaks = frameProxy.getNbrPeaks(scan_idx);
    //             if (numberOfPeaks == 0) {
    //                 continue;
    //             }

    //             // Create a new Spectrum object
    //             //TimsScan scan;
            
    //             // Populate Spectra for the current scan
    //             std::vector<double> mobilities;
    //             scanNumToOneOverK0(frame.getId(), {static_cast<double>(scan_idx)},mobilities);
    //             double mobility = mobilities.empty() ? 0.0 : mobilities[0];

    //             // Get the X (index) and Y (intensity) axis data
    //             auto x_axis = frameProxy.getScanX(scan_idx);
    //             auto y_axis = frameProxy.getScanY(scan_idx);

    //             // Convert X indices to m/z values
    //             std::vector<double> mz_values;
    //             std::vector<double> indices(x_axis.first, x_axis.second);
    //             indexToMz(frame.getId(), indices, mz_values);

    //             // Create and populate a Spectrum object
    //             TimsSpectrum spectrum(mobility);

    //             // size_t num_peaks = std::min(x_axis.size(), y_axis.size());
    //             for (size_t peak_idx = 0; peak_idx < numberOfPeaks; ++peak_idx) {
    //                 spectrum.addPeak(mz_values[peak_idx], y_axis.first[peak_idx]);
    //             }

    //             // Only add the Spectrum if it contains peaks
    //             if (spectrum.getNumberOfPeaks() > 0) {
    //                 //scan.addSpectrum(spectrum);
    //                 frame.addSpectrum(spectrum);
    //             }

    //             // Store the populated Scan object in the frame's vector
    //             //frame.addScan(scan);
    //         }
    //     }
    // }

    void TimsData::selectFrames(
        std::optional<int> min_id, std::optional<int> max_id,
        std::optional<double> min_time, std::optional<double> max_time,
        std::optional<std::string> polarity, std::optional<int> msms_type,
        std::optional<int> scan_mode, std::optional<int> tims_id,
        std::optional<double> min_max_intensity, std::optional<double> max_max_intensity,
        std::optional<double> min_summed_intensities, std::optional<double> max_summed_intensities,
        std::optional<int> min_num_peaks, std::optional<int> max_num_peaks) {

        // Check for invalid min/max value combinations and throw errors
        if (min_id && max_id && *min_id > *max_id) {
            throw std::invalid_argument("min_id cannot be greater than max_id");
        }
        if (min_time && max_time && *min_time > *max_time) {
            throw std::invalid_argument("min_time cannot be greater than max_time");
        }
        if (min_max_intensity && max_max_intensity && *min_max_intensity > *max_max_intensity) {
            throw std::invalid_argument("min_max_intensity cannot be greater than max_max_intensity");
        }
        if (min_summed_intensities && max_summed_intensities && *min_summed_intensities > *max_summed_intensities) {
            throw std::invalid_argument("min_summed_intensities cannot be greater than max_summed_intensities");
        }
        if (min_num_peaks && max_num_peaks && *min_num_peaks > *max_num_peaks) {
            throw std::invalid_argument("min_num_peaks cannot be greater than max_num_peaks");
        }
        
        std::stringstream filter;
        bool first_condition = true;

        auto append_condition = [&](const std::string& condition) {
            if (!first_condition) {
                filter << " AND ";
            }
            filter << condition;
            first_condition = false;
        };

        if (min_id) append_condition("Id >= " + std::to_string(*min_id));
        if (max_id) append_condition("Id <= " + std::to_string(*max_id));
        if (min_time) append_condition("Time >= " + std::to_string(*min_time));
        if (max_time) append_condition("Time <= " + std::to_string(*max_time));
        if (polarity) {
            if (*polarity == "+" || *polarity == "-") {
                append_condition("Polarity = '" + *polarity + "'");
            } else {
                throw std::invalid_argument("Polarity must be '+' or '-'");
            }
        }
        if (msms_type) append_condition("MsMsType = " + std::to_string(*msms_type));
        if (scan_mode) append_condition("ScanMode = " + std::to_string(*scan_mode));
        if (tims_id) append_condition("TimsId = " + std::to_string(*tims_id));
        if (min_max_intensity) append_condition("MaxIntensity >= " + std::to_string(*min_max_intensity));
        if (max_max_intensity) append_condition("MaxIntensity <= " + std::to_string(*max_max_intensity));
        if (min_summed_intensities) append_condition("SummedIntensities >= " + std::to_string(*min_summed_intensities));
        if (max_summed_intensities) append_condition("SummedIntensities <= " + std::to_string(*max_summed_intensities));
        if (min_num_peaks) append_condition("NumPeaks >= " + std::to_string(*min_num_peaks));
        if (max_num_peaks) append_condition("NumPeaks <= " + std::to_string(*max_num_peaks));

        // Convert the filter into a string
        std::string filter_str = filter.str();
    
        // Call populateFrames with the constructed filter
        populateFrames(filter_str);
    }

    void TimsData::populateSpectra(
        std::optional<double> min_mz,
        std::optional<double> max_mz,
        std::optional<double> min_mobility,
        std::optional<double> max_mobility) {
    
        // Ensure frames are populated first
        if (frames_.empty()) {
            throw std::runtime_error("Frames are not populated. Call populateFrames() first.");
        }

        // Loop through all frames
        for (auto& frame_pair : frames_) {
            Frame& frame = frame_pair.second;
            frame.clearSpectra();

            // Read scans for the current frame
            FrameProxy frameProxy = readScans(frame.getId(), 0, frame.getNumScans());

            for (uint32_t scan_idx = 0; scan_idx < frameProxy.getNbrScans(); ++scan_idx) {
                // Skip scans without peaks
                auto numberOfPeaks = frameProxy.getNbrPeaks(scan_idx);
                if (numberOfPeaks == 0) {
                    continue;
                }
        
                // Populate Spectra for the current scan
                std::vector<double> mobilities;
                scanNumToOneOverK0(frame.getId(), {static_cast<double>(scan_idx)}, mobilities);
                double mobility = mobilities.empty() ? 0.0 : mobilities[0];

                // Check mobility range
                if (min_mobility && mobility < *min_mobility) continue;
                if (max_mobility && mobility > *max_mobility) continue;

                // Get the X (index) and Y (intensity) axis data
                auto x_axis = frameProxy.getScanX(scan_idx);
                auto y_axis = frameProxy.getScanY(scan_idx);

                // Convert X indices to m/z values
                std::vector<double> mz_values;
                std::vector<double> indices(x_axis.first, x_axis.second);
                indexToMz(frame.getId(), indices, mz_values);

                // Create and populate a Spectrum object
                TimsSpectrum spectrum(mobility);

                for (size_t peak_idx = 0; peak_idx < numberOfPeaks; ++peak_idx) {
                    double mz = mz_values[peak_idx];
                    double intensity = y_axis.first[peak_idx];

                    // Check m/z range
                    if (min_mz && mz < *min_mz) continue;
                    if (max_mz && mz > *max_mz) continue;

                    spectrum.addPeak(mz, intensity);
                }

                // Only add the Spectrum if it contains peaks
                if (spectrum.getNumberOfPeaks() > 0) {
                    frame.addSpectrum(spectrum);
                }
            }
        }
    }

        
} // namespace timsdata


