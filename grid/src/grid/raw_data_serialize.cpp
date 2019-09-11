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
    Serialize::read_precursor_info(stream, &scan->precursor_information);
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
    Serialize::write_precursor_info(stream, scan.precursor_information);
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

bool IdentData::Serialize::read_spectrum_id(std::istream &stream,
                                            SpectrumId *spectrum_id) {
    Serialization::read_string(stream, &spectrum_id->id);
    Serialization::read_bool(stream, &spectrum_id->pass_threshold);
    Serialization::read_bool(stream, &spectrum_id->modifications);
    Serialization::read_string(stream, &spectrum_id->sequence);
    Serialization::read_string(stream, &spectrum_id->peptide_id);
    Serialization::read_uint64(stream, &spectrum_id->charge_state);
    Serialization::read_double(stream, &spectrum_id->theoretical_mz);
    Serialization::read_double(stream, &spectrum_id->experimental_mz);
    Serialization::read_double(stream, &spectrum_id->retention_time);
    Serialization::read_int64(stream, &spectrum_id->rank);
    return stream.good();
}

bool IdentData::Serialize::write_spectrum_id(std::ostream &stream,
                                             const SpectrumId &spectrum_id) {
    Serialization::write_string(stream, spectrum_id.id);
    Serialization::write_bool(stream, spectrum_id.pass_threshold);
    Serialization::write_bool(stream, spectrum_id.modifications);
    Serialization::write_string(stream, spectrum_id.sequence);
    Serialization::write_string(stream, spectrum_id.peptide_id);
    Serialization::write_uint64(stream, spectrum_id.charge_state);
    Serialization::write_double(stream, spectrum_id.theoretical_mz);
    Serialization::write_double(stream, spectrum_id.experimental_mz);
    Serialization::write_double(stream, spectrum_id.retention_time);
    Serialization::write_int64(stream, spectrum_id.rank);
    return stream.good();
}

bool IdentData::Serialize::read_cv_param(std::istream &stream,
                                         CVParam *cv_param) {
    Serialization::read_string(stream, &cv_param->name);
    Serialization::read_string(stream, &cv_param->accession);
    Serialization::read_string(stream, &cv_param->cv_ref);
    Serialization::read_string(stream, &cv_param->value);
    return stream.good();
}

bool IdentData::Serialize::write_cv_param(std::ostream &stream,
                                          const CVParam &cv_param) {
    Serialization::write_string(stream, cv_param.name);
    Serialization::write_string(stream, cv_param.accession);
    Serialization::write_string(stream, cv_param.cv_ref);
    Serialization::write_string(stream, cv_param.value);
    return stream.good();
}

bool IdentData::Serialize::read_peptide_mod(std::istream &stream,
                                            PeptideModification *peptide_mod) {
    Serialization::read_double(stream, &peptide_mod->monoisotopic_mass_delta);
    Serialization::read_double(stream, &peptide_mod->average_mass_delta);
    Serialization::read_string(stream, &peptide_mod->residues);
    Serialization::read_int64(stream, &peptide_mod->location);
    Serialization::read_vector<CVParam>(stream, &peptide_mod->cv_params,
                                        read_cv_param);
    return stream.good();
}

bool IdentData::Serialize::write_peptide_mod(
    std::ostream &stream, const PeptideModification &peptide_mod) {
    Serialization::write_double(stream, peptide_mod.monoisotopic_mass_delta);
    Serialization::write_double(stream, peptide_mod.average_mass_delta);
    Serialization::write_string(stream, peptide_mod.residues);
    Serialization::write_int64(stream, peptide_mod.location);
    Serialization::write_vector<CVParam>(stream, peptide_mod.cv_params,
                                         write_cv_param);
    return stream.good();
}

bool IdentData::Serialize::read_ident_data(std::istream &stream,
                                           IdentData *ident_data) {
    return stream.good();
}

bool IdentData::Serialize::write_ident_data(std::ostream &stream,
                                            const IdentData &ident_data) {
    return stream.good();
}
