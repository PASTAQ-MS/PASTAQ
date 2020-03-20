#ifndef RAWDATA_RAWDATASERIALIZE_HPP
#define RAWDATA_RAWDATASERIALIZE_HPP

#include <iostream>

#include "raw_data/raw_data.hpp"

// This namespace groups the functions used to serialize RawData data structures
// into a binary stream.
namespace RawData::Serialize {

// RawData::Scan
bool read_scan(std::istream &stream, Scan *scan);
bool write_scan(std::ostream &stream, const Scan &scan);

// RawData::PrecursorInformation
bool read_precursor_info(std::istream &stream,
                         PrecursorInformation *precursor_info);
bool write_precursor_info(std::ostream &stream,
                          const PrecursorInformation &precursor_info);

// RawData::RawData
bool read_raw_data(std::istream &stream, RawData *raw_data);
bool write_raw_data(std::ostream &stream, const RawData &raw_data);

}  // namespace RawData::Serialize

// This namespace groups the functions used to serialize IdentData data
// structures into a binary stream.
namespace IdentData::Serialize {

// IdentData::SpectrumMatch
bool read_spectrum_match(std::istream &stream, SpectrumMatch *spectrum_match);
bool write_spectrum_match(std::ostream &stream,
                          const SpectrumMatch &spectrum_match);

// IdentData::DBSequence
bool read_db_sequence(std::istream &stream, DBSequence *db_sequence);
bool write_db_sequence(std::ostream &stream, const DBSequence &db_sequence);

// IdentData::PeptideModification
bool read_peptide_mod(std::istream &stream, PeptideModification *peptide_mod);
bool write_peptide_mod(std::ostream &stream,
                       const PeptideModification &peptide_mod);

// IdentData::Peptide
bool read_peptide(std::istream &stream, Peptide *peptide);
bool write_peptide(std::ostream &stream, const Peptide &peptide);

// IdentData::PeptideEvidence
bool read_peptide_evidence(std::istream &stream,
                           PeptideEvidence *peptide_evidence);
bool write_peptide_evidence(std::ostream &stream,
                            const PeptideEvidence &peptide_evidence);

// IdentData::IdentData
bool read_ident_data(std::istream &stream, IdentData *ident_data);
bool write_ident_data(std::ostream &stream, const IdentData &ident_data);
}  // namespace IdentData::Serialize

#endif /* RAWDATA_RAWDATASERIALIZE_HPP */
