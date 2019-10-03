#include "centroid/centroid_serialize.hpp"
#include "utils/serialization.hpp"

bool Centroid::Serialize::read_peak(std::istream &stream,
                                    Centroid::Peak *peak) {
    Serialization::read_uint64(stream, &peak->id);
    Serialization::read_double(stream, &peak->local_max_mz);
    Serialization::read_double(stream, &peak->local_max_rt);
    Serialization::read_double(stream, &peak->local_max_height);
    Serialization::read_double(stream, &peak->roi_min_mz);
    Serialization::read_double(stream, &peak->roi_max_mz);
    Serialization::read_double(stream, &peak->roi_min_rt);
    Serialization::read_double(stream, &peak->roi_max_rt);
    Serialization::read_double(stream, &peak->raw_roi_mean_mz);
    Serialization::read_double(stream, &peak->raw_roi_mean_rt);
    Serialization::read_double(stream, &peak->raw_roi_sigma_mz);
    Serialization::read_double(stream, &peak->raw_roi_sigma_rt);
    Serialization::read_double(stream, &peak->raw_roi_skewness_mz);
    Serialization::read_double(stream, &peak->raw_roi_skewness_rt);
    Serialization::read_double(stream, &peak->raw_roi_kurtosis_mz);
    Serialization::read_double(stream, &peak->raw_roi_kurtosis_rt);
    Serialization::read_double(stream, &peak->raw_roi_max_height);
    Serialization::read_double(stream, &peak->raw_roi_total_intensity);
    Serialization::read_uint64(stream, &peak->raw_roi_num_points);
    Serialization::read_uint64(stream, &peak->raw_roi_num_scans);
    return stream.good();
}

bool Centroid::Serialize::write_peak(std::ostream &stream,
                                     const Centroid::Peak &peak) {
    Serialization::write_uint64(stream, peak.id);
    Serialization::write_double(stream, peak.local_max_mz);
    Serialization::write_double(stream, peak.local_max_rt);
    Serialization::write_double(stream, peak.local_max_height);
    Serialization::write_double(stream, peak.roi_min_mz);
    Serialization::write_double(stream, peak.roi_max_mz);
    Serialization::write_double(stream, peak.roi_min_rt);
    Serialization::write_double(stream, peak.roi_max_rt);
    Serialization::write_double(stream, peak.raw_roi_mean_mz);
    Serialization::write_double(stream, peak.raw_roi_mean_rt);
    Serialization::write_double(stream, peak.raw_roi_sigma_mz);
    Serialization::write_double(stream, peak.raw_roi_sigma_rt);
    Serialization::write_double(stream, peak.raw_roi_skewness_mz);
    Serialization::write_double(stream, peak.raw_roi_skewness_rt);
    Serialization::write_double(stream, peak.raw_roi_kurtosis_mz);
    Serialization::write_double(stream, peak.raw_roi_kurtosis_rt);
    Serialization::write_double(stream, peak.raw_roi_max_height);
    Serialization::write_double(stream, peak.raw_roi_total_intensity);
    Serialization::write_uint64(stream, peak.raw_roi_num_points);
    Serialization::write_uint64(stream, peak.raw_roi_num_scans);
    return stream.good();
}

bool Centroid::Serialize::read_peaks(std::istream &stream,
                                       std::vector<Centroid::Peak> *peaks) {
    uint64_t num_peaks = 0;
    Serialization::read_uint64(stream, &num_peaks);
    peaks->resize(num_peaks);
    for (auto &peak : *peaks) {
        if (!Centroid::Serialize::read_peak(stream, &peak)) {
            return false;
        }
    }
    return stream.good();
}

bool Centroid::Serialize::write_peaks(
    std::ostream &stream, const std::vector<Centroid::Peak> &peaks) {
    Serialization::write_uint64(stream, peaks.size());
    for (const auto &peak : peaks) {
        if (!Centroid::Serialize::write_peak(stream, peak)) {
            return false;
        }
    }
    return stream.good();
}
