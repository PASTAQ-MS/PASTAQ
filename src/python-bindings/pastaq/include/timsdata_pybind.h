#ifndef TIMSDATACPP_PYBIND_H
#define TIMSDATACPP_PYBIND_H

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>  // For handling std::unique_ptr and arrays
#include "timsdatacpp.h"  // Include core C++ class

namespace py = pybind11;

void register_tims_data(py::module_ &m);

class PyFrameProxy {
public:
    PyFrameProxy(uint32_t num_scans, std::shared_ptr<uint32_t[]> data)
        : num_scans(num_scans), data(std::move(data)), scan_offsets(num_scans + 1) {
        scan_offsets[0] = 0;
        std::partial_sum(this->data.get(), this->data.get() + num_scans, scan_offsets.begin() + 1);
    }

    size_t getNbrScans() const {
        return num_scans;
    }

    size_t getTotalNbrPeaks() const {
        return scan_offsets.back();
    }

    size_t getNbrPeaks(size_t scan_num) const {
        if (scan_num >= num_scans)
            throw std::out_of_range("Invalid scan number");
        return scan_offsets[scan_num + 1] - scan_offsets[scan_num];
    }

    std::vector<uint32_t> get_data() const {
        return std::vector<uint32_t>(data.get(), data.get() + num_scans);
    }

    std::vector<uint32_t> getScanX(size_t scan_num) const {
        auto range = makeRange(scan_num, 0);
        return std::vector<uint32_t>(range.first, range.second);
    }

    std::vector<uint32_t> getScanY(size_t scan_num) const {
        auto range = makeRange(scan_num, data[scan_num]);
        return std::vector<uint32_t>(range.first, range.second);
    }
private:
    uint32_t num_scans;
    std::shared_ptr<uint32_t[]> data;

    std::vector<uint32_t> scan_offsets;

    void throwIfInvalidScanNumber(size_t scan_num) const {
        if (scan_num >= num_scans)
            throw std::out_of_range("Invalid scan number");
    }

    std::pair<const uint32_t*, const uint32_t*> makeRange(size_t scan_num, size_t offset) const {
        throwIfInvalidScanNumber(scan_num);
        const uint32_t* start = data.get() + num_scans + 2 * scan_offsets[scan_num] + offset;
        const uint32_t* end = start + data[scan_num]; // Adjust end pointer like the original function

        return {start, end};
    }
};

#endif // TIMSDATACPP_PYBIND_H