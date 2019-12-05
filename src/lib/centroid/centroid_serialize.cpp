#include "centroid/centroid_serialize.hpp"
#include "utils/serialization.hpp"

bool Centroid::Serialize::read_peak(std::istream &stream,
                                    Centroid::Peak *peak) {
    Serialization::read_uint64(stream, &peak->id);
    Serialization::read_double(stream, &peak->local_max_mz);
    Serialization::read_double(stream, &peak->local_max_rt);
    Serialization::read_double(stream, &peak->local_max_height);
    Serialization::read_double(stream, &peak->rt_delta);
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
    Serialization::read_double(stream, &peak->fitted_height);
    Serialization::read_double(stream, &peak->fitted_mz);
    Serialization::read_double(stream, &peak->fitted_rt);
    Serialization::read_double(stream, &peak->fitted_sigma_mz);
    Serialization::read_double(stream, &peak->fitted_sigma_rt);
    Serialization::read_double(stream, &peak->fitted_volume);
    return stream.good();
}

bool Centroid::Serialize::write_peak(std::ostream &stream,
                                     const Centroid::Peak &peak) {
    Serialization::write_uint64(stream, peak.id);
    Serialization::write_double(stream, peak.local_max_mz);
    Serialization::write_double(stream, peak.local_max_rt);
    Serialization::write_double(stream, peak.local_max_height);
    Serialization::write_double(stream, peak.rt_delta);
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
    Serialization::write_double(stream, peak.fitted_height);
    Serialization::write_double(stream, peak.fitted_mz);
    Serialization::write_double(stream, peak.fitted_rt);
    Serialization::write_double(stream, peak.fitted_sigma_mz);
    Serialization::write_double(stream, peak.fitted_sigma_rt);
    Serialization::write_double(stream, peak.fitted_volume);
    return stream.good();
}

bool Centroid::Serialize::read_peaks(std::istream &stream,
                                     std::vector<Centroid::Peak> *peaks) {
    return Serialization::read_vector<Centroid::Peak>(
        stream, peaks, Centroid::Serialize::read_peak);
}

bool Centroid::Serialize::write_peaks(
    std::ostream &stream, const std::vector<Centroid::Peak> &peaks) {
    return Serialization::write_vector<Centroid::Peak>(
        stream, peaks, Centroid::Serialize::write_peak);
}
