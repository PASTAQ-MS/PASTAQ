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
#include <utility>
#include <optional>

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

    // Forward Declaration
    class FrameProxy;
    class Spectrum;
    class Scan;
    class Frame;
    class TimsData;

    class TimsSpectrum {
        public:
            TimsSpectrum(double mobility) : mobility_(mobility) {}

            void addPeak(double mzValue, double intensityValue) {
                mz_.push_back(mzValue);
                intensity_.push_back(intensityValue);
            }

            const std::vector<double>& getMz() const {
                return mz_;
            }

            const std::vector<double>& getIntensity() const {
                return intensity_;
            }

            double getMobility() const {
                return mobility_;
            }

            size_t getNumberOfPeaks() const {
                return mz_.size();
            }

        private:
            double mobility_;
            std::vector<double> mz_;
            std::vector<double> intensity_;
    };

    // Scan Class
    // class TimsScan {
    //     public:
    //         void addSpectrum(const TimsSpectrum& spectrum) {
    //             spectra_.push_back(spectrum);
    //         }

    //         const std::vector<TimsSpectrum>& getSpectra() const {
    //             return spectra_;
    //         }

    //         size_t getNumberOfSpectra() const {
    //             return spectra_.size();
    //         }

    //     private:
    //         std::vector<TimsSpectrum> spectra_;
    // };

    // Frame Class
    class Frame {
        public:
            void addSpectrum(const TimsSpectrum& spectrum) {
                spectra_.push_back(spectrum);
            }

            void clearSpectra() {
                spectra_.clear();
            }

            const std::vector<TimsSpectrum>& getSpectra() const {
                return spectra_;
            }

            size_t getNumberOfScans() const {
                return spectra_.size();
            }

            // Setters 
            void setId(int64_t id) { Id = id; }
            void setTime(double time) { Time = time; }
            void setPolarity(const std::string& polarity) { Polarity = polarity; }
            void setScanMode(int scanMode) { ScanMode = scanMode; }
            void setMsMsType(int msMsType) { MsMsType = msMsType; }
            void setTimsId(int32_t timsId) { TimsId = timsId; }
            void setMaxIntensity(uint32_t maxIntensity) { MaxIntensity = maxIntensity; }
            void setSummedIntensities(uint64_t summedIntensities) { SummedIntensities = summedIntensities; }
            void setNumScans(uint32_t numScans) { NumScans = numScans; }
            void setNumPeaks(uint32_t numPeaks) { NumPeaks = numPeaks; }
            void setMzCalibration(double mzCalibration) { MzCalibration = mzCalibration; }
            void setT1(double t1) { T1 = t1; }
            void setT2(double t2) { T2 = t2; }
            void setTimsCalibration(double timsCalibration) { TimsCalibration = timsCalibration; }
            void setPropertyGroup(int propertyGroup) { PropertyGroup = propertyGroup; }
            void setAccumulationTime(double accumulationTime) { AccumulationTime = accumulationTime; }
            void setRampTime(double rampTime) { RampTime = rampTime; }
            void setPressure(double pressure) { Pressure = pressure; }

            // Getters
            int64_t getId() const { return Id; }
            double getTime() const { return Time; }
            const std::string& getPolarity() const { return Polarity; }
            int getScanMode() const { return ScanMode; }
            int getMsMsType() const { return MsMsType; }
            int32_t getTimsId() const { return TimsId; }
            uint32_t getMaxIntensity() const { return MaxIntensity; }
            uint64_t getSummedIntensities() const { return SummedIntensities; }
            uint32_t getNumScans() const { return NumScans; }
            uint32_t getNumPeaks() const { return NumPeaks; }
            double getMzCalibration() const { return MzCalibration; }
            double getT1() const { return T1; }
            double getT2() const { return T2; }
            double getTimsCalibration() const { return TimsCalibration; }
            int getPropertyGroup() const { return PropertyGroup; }
            double getAccumulationTime() const { return AccumulationTime; }
            double getRampTime() const { return RampTime; }
            double getPressure() const { return Pressure; }

        private:
            // std::vector<TimsScan> scans_;
            std::vector<TimsSpectrum> spectra_;

            int64_t Id;
            double Time;
            std::string Polarity;
            int ScanMode;
            int MsMsType;
            int32_t TimsId;
            uint32_t MaxIntensity;
            uint64_t SummedIntensities;
            uint32_t NumScans;
            uint32_t NumPeaks;
            double MzCalibration;
            double T1;
            double T2;
            double TimsCalibration;
            int PropertyGroup;
            double AccumulationTime;
            double RampTime;
            double Pressure;
    };

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
        void populateFrames(const std::string& filter = "");
        const std::map<int64_t, Frame>& getFrames() const {
            return frames_;
        }

        void selectFrames(
            std::optional<int> min_id, std::optional<int> max_id,
            std::optional<double> min_time, std::optional<double> max_time,
            std::optional<std::string> polarity, std::optional<int> msms_type,
            std::optional<int> scan_mode, std::optional<int> tims_id,
            std::optional<double> min_max_intensity, std::optional<double> max_max_intensity,
            std::optional<double> min_summed_intensities, std::optional<double> max_summed_intensities,
            std::optional<int> min_num_peaks, std::optional<int> max_num_peaks);

        void populateSpectra(std::optional<double> min_mz = std::nullopt,
            std::optional<double> max_mz = std::nullopt,
            std::optional<double> min_mobility = std::nullopt,
            std::optional<double> max_mobility = std::nullopt);

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

        void doTransformation(int64_t frame_id, 
                              const std::vector<double>& in,
                              std::vector<double>& out, 
                              BdalTimsConversionFunction* func);
        std::string tdfFile;             // Path to the .tdf SQLite file
        mutable CppSQLite3DB db;         // Database connection object

        // Cache for frames
        std::map<int64_t, Frame> frames_;
    };

} // namespace timsdata

#endif // DE_BDAL_TIMSDATA_CPP_H
