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
} // namespace timsdata
