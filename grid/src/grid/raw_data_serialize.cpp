#include "raw_data_serialize.hpp"
#include "utils/serialization.hpp"

bool RawData::Serialize::read_precursor_info(
    std::istream &stream, PrecursorInformation *precursor_info) {
    Serialization::read_uint64(stream, &precursor_info->scan_number);
    Serialization::read_uint8(stream, &precursor_info->charge);
    Serialization::read_double(stream, &precursor_info->mz);
    Serialization::read_double(stream, &precursor_info->intensity);
    uint8_t activation_method = 0;
    Serialization::read_uint8(stream, &activation_method);
    precursor_info->activation_method =
        static_cast<ActivationMethod>(activation_method);
    Serialization::read_double(stream, &precursor_info->window_wideness);
    return stream.good();
}

bool RawData::Serialize::write_precursor_info(
    std::ostream &stream, const PrecursorInformation &precursor_info) {
    Serialization::write_uint64(stream, precursor_info.scan_number);
    Serialization::write_uint8(stream, precursor_info.charge);
    Serialization::write_double(stream, precursor_info.mz);
    Serialization::write_double(stream, precursor_info.intensity);
    Serialization::write_uint8(stream, precursor_info.activation_method);
    Serialization::write_double(stream, precursor_info.window_wideness);
    return stream.good();
}

bool RawData::Serialize::read_scan(std::istream &stream, Scan *scan) {
    Serialization::read_uint64(stream, &scan->scan_number);
    Serialization::read_uint64(stream, &scan->ms_level);
    Serialization::read_uint64(stream, &scan->num_points);
    Serialization::read_double(stream, &scan->retention_time);
    scan->mz = std::vector<double>(scan->num_points);
    scan->intensity = std::vector<double>(scan->num_points);
    for (size_t i = 0; i < scan->num_points; ++i) {
        Serialization::read_double(stream, &scan->mz[i]);
        Serialization::read_double(stream, &scan->intensity[i]);
    }
    // TODO: Precursor information if pointer is present.
    return stream.good();
}

bool RawData::Serialize::write_scan(std::ostream &stream, const Scan &scan) {
    Serialization::write_uint64(stream, scan.scan_number);
    Serialization::write_uint64(stream, scan.ms_level);
    Serialization::write_uint64(stream, scan.num_points);
    Serialization::write_double(stream, scan.retention_time);
    for (size_t i = 0; i < scan.num_points; ++i) {
        Serialization::write_double(stream, scan.mz[i]);
        Serialization::write_double(stream, scan.intensity[i]);
    }
    // TODO: Precursor information if pointer is present.
    return stream.good();
}

bool RawData::Serialize::read_raw_data(std::istream &stream,
                                       RawData *raw_data) {
    uint8_t instrument_type = 0;
    Serialization::read_uint8(stream, &instrument_type);
    raw_data->instrument_type = static_cast<Instrument::Type>(instrument_type);
    Serialization::read_double(stream, &raw_data->min_mz);
    Serialization::read_double(stream, &raw_data->max_mz);
    Serialization::read_double(stream, &raw_data->min_rt);
    Serialization::read_double(stream, &raw_data->max_rt);
    Serialization::read_double(stream, &raw_data->resolution_ms1);
    Serialization::read_double(stream, &raw_data->resolution_msn);
    Serialization::read_double(stream, &raw_data->reference_mz);
    Serialization::read_double(stream, &raw_data->fwhm_rt);
    uint64_t num_scans = 0;
    Serialization::read_uint64(stream, &num_scans);
    raw_data->scans = std::vector<Scan>(num_scans);
    raw_data->retention_times = std::vector<double>(num_scans);
    raw_data->total_ion_chromatogram = std::vector<double>(num_scans);
    raw_data->base_peak_chromatogram = std::vector<double>(num_scans);
    for (size_t i = 0; i < num_scans; ++i) {
        Serialize::read_scan(stream, &raw_data->scans[i]);
        Serialization::read_double(stream, &raw_data->retention_times[i]);
        Serialization::read_double(stream,
                                   &raw_data->total_ion_chromatogram[i]);
        Serialization::read_double(stream,
                                   &raw_data->base_peak_chromatogram[i]);
    }
    // TODO: Precursor information if pointer is present.
    return stream.good();
}

bool RawData::Serialize::write_raw_data(std::ostream &stream,
                                        const RawData &raw_data) {
    Serialization::write_uint8(stream, raw_data.instrument_type);
    Serialization::write_double(stream, raw_data.min_mz);
    Serialization::write_double(stream, raw_data.max_mz);
    Serialization::write_double(stream, raw_data.min_rt);
    Serialization::write_double(stream, raw_data.max_rt);
    Serialization::write_double(stream, raw_data.resolution_ms1);
    Serialization::write_double(stream, raw_data.resolution_msn);
    Serialization::write_double(stream, raw_data.reference_mz);
    Serialization::write_double(stream, raw_data.fwhm_rt);
    uint64_t num_scans = raw_data.scans.size();
    Serialization::write_uint64(stream, num_scans);
    for (size_t i = 0; i < num_scans; ++i) {
        Serialize::write_scan(stream, raw_data.scans[i]);
        Serialization::write_double(stream, raw_data.retention_times[i]);
        Serialization::write_double(stream, raw_data.total_ion_chromatogram[i]);
        Serialization::write_double(stream, raw_data.base_peak_chromatogram[i]);
    }
    return stream.good();
}