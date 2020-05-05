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
        static_cast<ActivationMethod::Type>(activation_method);
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
    uint8_t polarity = polarity;
    Serialization::read_uint8(stream, &polarity);
    scan->polarity = static_cast<Polarity::Type>(polarity);
    Serialization::read_double(stream, &scan->max_intensity);
    Serialization::read_double(stream, &scan->total_intensity);
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
    Serialization::write_uint8(stream, scan.polarity);
    Serialization::write_double(stream, scan.max_intensity);
    Serialization::write_double(stream, scan.total_intensity);
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
    for (size_t i = 0; i < num_scans; ++i) {
        Serialize::read_scan(stream, &raw_data->scans[i]);
        Serialization::read_double(stream, &raw_data->retention_times[i]);
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
    }
    return stream.good();
}

bool IdentData::Serialize::read_spectrum_match(std::istream &stream,
                                               SpectrumMatch *spectrum_match) {
    Serialization::read_string(stream, &spectrum_match->id);
    Serialization::read_bool(stream, &spectrum_match->pass_threshold);
    Serialization::read_string(stream, &spectrum_match->match_id);
    Serialization::read_uint8(stream, &spectrum_match->charge_state);
    Serialization::read_double(stream, &spectrum_match->theoretical_mz);
    Serialization::read_double(stream, &spectrum_match->experimental_mz);
    Serialization::read_double(stream, &spectrum_match->retention_time);
    Serialization::read_uint64(stream, &spectrum_match->rank);
    return stream.good();
}

bool IdentData::Serialize::write_spectrum_match(
    std::ostream &stream, const SpectrumMatch &spectrum_match) {
    Serialization::write_string(stream, spectrum_match.id);
    Serialization::write_bool(stream, spectrum_match.pass_threshold);
    Serialization::write_string(stream, spectrum_match.match_id);
    Serialization::write_uint8(stream, spectrum_match.charge_state);
    Serialization::write_double(stream, spectrum_match.theoretical_mz);
    Serialization::write_double(stream, spectrum_match.experimental_mz);
    Serialization::write_double(stream, spectrum_match.retention_time);
    Serialization::write_uint64(stream, spectrum_match.rank);
    return stream.good();
}

bool IdentData::Serialize::read_db_sequence(std::istream &stream,
                                            DBSequence *db_sequence) {
    Serialization::read_string(stream, &db_sequence->id);
    Serialization::read_string(stream, &db_sequence->accession);
    Serialization::read_string(stream, &db_sequence->db_reference);
    Serialization::read_string(stream, &db_sequence->description);
    return stream.good();
}

bool IdentData::Serialize::write_db_sequence(std::ostream &stream,
                                             const DBSequence &db_sequence) {
    Serialization::write_string(stream, db_sequence.id);
    Serialization::write_string(stream, db_sequence.accession);
    Serialization::write_string(stream, db_sequence.db_reference);
    Serialization::write_string(stream, db_sequence.description);
    return stream.good();
}

bool IdentData::Serialize::read_peptide_mod(std::istream &stream,
                                            PeptideModification *peptide_mod) {
    Serialization::read_double(stream, &peptide_mod->monoisotopic_mass_delta);
    Serialization::read_double(stream, &peptide_mod->average_mass_delta);
    Serialization::read_string(stream, &peptide_mod->residues);
    Serialization::read_int64(stream, &peptide_mod->location);
    Serialization::read_vector<std::string>(stream, &peptide_mod->id,
                                            Serialization::read_string);
    return stream.good();
}

bool IdentData::Serialize::write_peptide_mod(
    std::ostream &stream, const PeptideModification &peptide_mod) {
    Serialization::write_double(stream, peptide_mod.monoisotopic_mass_delta);
    Serialization::write_double(stream, peptide_mod.average_mass_delta);
    Serialization::write_string(stream, peptide_mod.residues);
    Serialization::write_int64(stream, peptide_mod.location);
    Serialization::write_vector<std::string>(stream, peptide_mod.id,
                                             Serialization::write_string);
    return stream.good();
}

bool IdentData::Serialize::read_peptide(std::istream &stream,
                                        Peptide *peptide) {
    Serialization::read_string(stream, &peptide->id);
    Serialization::read_string(stream, &peptide->sequence);
    Serialization::read_vector<PeptideModification>(
        stream, &peptide->modifications, read_peptide_mod);
    return stream.good();
}

bool IdentData::Serialize::write_peptide(std::ostream &stream,
                                         const Peptide &peptide) {
    Serialization::write_string(stream, peptide.id);
    Serialization::write_string(stream, peptide.sequence);
    Serialization::write_vector<PeptideModification>(
        stream, peptide.modifications, write_peptide_mod);
    return stream.good();
}

bool IdentData::Serialize::read_peptide_evidence(
    std::istream &stream, PeptideEvidence *peptide_evidence) {
    Serialization::read_string(stream, &peptide_evidence->id);
    Serialization::read_string(stream, &peptide_evidence->db_sequence_id);
    Serialization::read_string(stream, &peptide_evidence->peptide_id);
    Serialization::read_bool(stream, &peptide_evidence->decoy);
    return stream.good();
}

bool IdentData::Serialize::write_peptide_evidence(
    std::ostream &stream, const PeptideEvidence &peptide_evidence) {
    Serialization::write_string(stream, peptide_evidence.id);
    Serialization::write_string(stream, peptide_evidence.db_sequence_id);
    Serialization::write_string(stream, peptide_evidence.peptide_id);
    Serialization::write_bool(stream, peptide_evidence.decoy);
    return stream.good();
}

bool IdentData::Serialize::read_ident_data(std::istream &stream,
                                           IdentData *ident_data) {
    Serialization::read_vector<DBSequence>(stream, &ident_data->db_sequences,
                                           read_db_sequence);
    Serialization::read_vector<Peptide>(stream, &ident_data->peptides,
                                        read_peptide);
    Serialization::read_vector<SpectrumMatch>(
        stream, &ident_data->spectrum_matches, read_spectrum_match);
    Serialization::read_vector<PeptideEvidence>(
        stream, &ident_data->peptide_evidence, read_peptide_evidence);
    return stream.good();
}

bool IdentData::Serialize::write_ident_data(std::ostream &stream,
                                            const IdentData &ident_data) {
    Serialization::write_vector<DBSequence>(stream, ident_data.db_sequences,
                                            write_db_sequence);
    Serialization::write_vector<Peptide>(stream, ident_data.peptides,
                                         write_peptide);
    Serialization::write_vector<SpectrumMatch>(
        stream, ident_data.spectrum_matches, write_spectrum_match);
    Serialization::write_vector<PeptideEvidence>(
        stream, ident_data.peptide_evidence, write_peptide_evidence);
    return stream.good();
}
