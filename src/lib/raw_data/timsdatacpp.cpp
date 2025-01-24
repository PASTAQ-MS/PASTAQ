#include <stdexcept>
#include <string>
#include <cstdint>
#include <numeric>
#include <vector>
#include <limits>
#include <memory>
#include <utility> // For std::move
#include <iostream>

// #include <cstddef>
#include "CppSQLite3.h"

#include "timsdata.h"
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

} // namespace timsdata
