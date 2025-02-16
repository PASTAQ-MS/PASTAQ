#ifndef DE_BDAL_TIMSDATA_CPP_H
#define DE_BDAL_TIMSDATA_CPP_H

#include <map>
#include <sstream>
#include <cstddef>     // For size_t
#include <stdexcept>
#include <string>
#include <cstdint>
#include <numeric>
#include <vector>
#include <limits>
#include <memory>
#include <variant>

#include "CppSQLite3.h"
#include "timsdata.h" // fundamental C API

        #define BDAL_TIMS_DEFINE_CONVERSION_FUNCTION(CPPNAME, CNAME) \
        void CPPNAME ( \
            int64_t frame_id,               /**< frame index */ \
            const std::vector<double> & in, /**< vector of input values (can be empty) */ \
            std::vector<double> & out )     /**< vector of corresponding output values (will be resized automatically) */ \
        { \
            doTransformation(frame_id, in, out, CNAME); \
        }


namespace timsdata
{

    class FrameProxy
    {
    public:
        using FrameIteratorRange = std::pair<const uint32_t*, const uint32_t*>;

        FrameProxy(size_t num_scans_, std::shared_ptr<uint32_t[]> pData_)
        : num_scans(num_scans_)
        , pData(std::move(pData_))
        , scan_offsets(num_scans_ + 1)
        {
            scan_offsets[0] = 0;
            std::partial_sum(pData.get(), pData.get() + num_scans, scan_offsets.begin() + 1);
        }

        size_t getNbrScans() const;
        size_t getTotalNbrPeaks() const;
        size_t getNbrPeaks(size_t scan_num) const;
        FrameIteratorRange getScanX(size_t scan_num) const;
        FrameIteratorRange getScanY(size_t scan_num) const;

        const std::shared_ptr<uint32_t[]>& getData() const {
            return pData;
        }

    private:
        const size_t num_scans;
        const std::shared_ptr<uint32_t[]> pData;
        std::vector<uint32_t> scan_offsets;

        void throwIfInvalidScanNumber(size_t scan_num) const;
        FrameIteratorRange makeRange(size_t scan_num, size_t offset) const;
    };

    std::string getLastError();
    void throwLastError();

    class TimsData
    {
    public:
        explicit TimsData(const std::string& analysis_directory_name, bool use_recalibration = false, 
                            pressure_compensation_strategy pressure_compensation = AnalyisGlobalPressureCompensation);
        ~TimsData();

        TimsData(const TimsData&) = delete;
        TimsData operator=(const TimsData&) = delete;

        TimsData(TimsData&& other) = default;
        TimsData& operator=(TimsData&&) = default;

        uint64_t getHandle() const;

        FrameProxy readScans(int64_t frame_id, uint32_t scan_begin, uint32_t scan_end);

        uint32_t getNumberOfFrames() const;
        std::pair<std::vector<std::string>, std::vector<std::vector<std::variant<int64_t, double, std::string>>>> getFramesTable() const;
        std::string getTdfFile() const {
            return tdfFile;
        }

        BDAL_TIMS_DEFINE_CONVERSION_FUNCTION(indexToMz, tims_index_to_mz)
        BDAL_TIMS_DEFINE_CONVERSION_FUNCTION(mzToIndex, tims_mz_to_index)
        BDAL_TIMS_DEFINE_CONVERSION_FUNCTION(scanNumToOneOverK0, tims_scannum_to_oneoverk0)
        BDAL_TIMS_DEFINE_CONVERSION_FUNCTION(oneOverK0ToScanNum, tims_oneoverk0_to_scannum)
        BDAL_TIMS_DEFINE_CONVERSION_FUNCTION(scanNumToVoltage, tims_scannum_to_voltage)
        BDAL_TIMS_DEFINE_CONVERSION_FUNCTION(voltageToScanNum, tims_voltage_to_scannum)

        std::vector<std::string> get_tables() const;
        std::string get_schema(const std::string& table_name) const;
        std::vector<std::map<std::string, std::variant<int64_t, double, std::string>>> query(const std::string& sql_query) const;
        void relationships() const;

    private:
        uint64_t handle;
        size_t initial_frame_buffer_size;

        void doTransformation(int64_t frame_id, const std::vector<double>& in, std::vector<double>& out, 
                                BdalTimsConversionFunction* func);
        std::string tdfFile;             // Path to the .tdf SQLite file
        mutable CppSQLite3DB db;         // Database connection object
    };

} // namespace timsdata

#endif // DE_BDAL_TIMSDATA_CPP_H
