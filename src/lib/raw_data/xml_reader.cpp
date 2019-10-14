#include <regex>
#include <sstream>

#include "utils/base64.hpp"
#include "xml_reader.hpp"

std::optional<RawData::RawData> XmlReader::read_mzxml(
    std::istream &stream, double min_mz, double max_mz, double min_rt,
    double max_rt, Instrument::Type instrument_type, double resolution_ms1,
    double resolution_msn, double reference_mz, Polarity::Type polarity,
    size_t ms_level) {
    RawData::RawData raw_data = {};
    raw_data.instrument_type = instrument_type;
    raw_data.min_mz = std::numeric_limits<double>::infinity();
    raw_data.max_mz = -std::numeric_limits<double>::infinity();
    raw_data.min_rt = std::numeric_limits<double>::infinity();
    raw_data.max_rt = -std::numeric_limits<double>::infinity();
    raw_data.resolution_ms1 = resolution_ms1;
    raw_data.resolution_msn = resolution_msn;
    raw_data.reference_mz = reference_mz;
    raw_data.fwhm_rt = 0;  // TODO(alex): Should this be passed as well?
    raw_data.scans = {};
    raw_data.retention_times = {};
    // TODO(alex): Can we automatically detect the instrument type and set
    // resolution from the header?
    size_t previous_ms1_scan_number = 0;
    while (stream.good()) {
        auto tag = XmlReader::read_tag(stream);
        if (!tag) {
            continue;
        }
        if (tag.value().name == "scan" && !tag.value().closed) {
            RawData::Scan scan;

            auto scan_attributes = tag.value().attributes;

            // Find scan number.
            if (scan_attributes.find("num") == scan_attributes.end()) {
                return std::nullopt;
            }
            scan.scan_number = std::stoi(scan_attributes["num"]);

            // Find polarity.
            if (scan_attributes.find("polarity") != scan_attributes.end()) {
                if (scan_attributes["polarity"] == "+") {
                    scan.polarity = Polarity::POSITIVE;
                } else if (scan_attributes["polarity"] == "-") {
                    scan.polarity = Polarity::NEGATIVE;
                } else {
                    scan.polarity = Polarity::BOTH;
                }
                if (polarity != Polarity::BOTH && scan.polarity != polarity) {
                    continue;
                }
            }

            // Find MS level.
            if (scan_attributes.find("msLevel") == scan_attributes.end()) {
                return std::nullopt;
            }
            size_t scan_ms_level = std::stoi(scan_attributes["msLevel"]);
            if (scan_ms_level == 1) {
                previous_ms1_scan_number = scan.scan_number;
            }
            if (scan_ms_level != ms_level) {
                continue;
            }
            scan.ms_level = scan_ms_level;

            // Find the number of m/z-intensity pairs in the scan.
            if (scan_attributes.find("peaksCount") == scan_attributes.end()) {
                return std::nullopt;
            }
            size_t num_points = std::stoi(scan_attributes["peaksCount"]);

            // Extract the retention time.
            if (scan_attributes.find("retentionTime") ==
                scan_attributes.end()) {
                // NOTE(alex): On the spec, the retention time attribute is
                // optional, however, we do require it.
                return std::nullopt;
            }

            // The time is in xs:duration units. Here we are only accounting
            // for the data as stored in seconds, minutes and hours. It is
            // unlikely that we are going to need to parse the days, months
            // and years. For more information about the format see:
            //    https://www.ibm.com/support/knowledgecenter/en/ssw_ibm_i_72/rzasp/rzasp_xsduration.htm
            std::regex rt_regex(
                R"(P.*T(?:([[:digit:]]+)H)?(?:([[:digit:]]+)M)?(?:([[:digit:]]+\.?[[:digit:]]*)S))");
            std::smatch matches;
            if (!std::regex_search(scan_attributes["retentionTime"], matches,
                                   rt_regex) ||
                matches.size() != 4) {
                return std::nullopt;
            }
            double retention_time = std::stod(matches[3]);
            if (matches[2] != "") {
                retention_time += std::stod(matches[2]) * 60;
            }
            if (matches[1] != "") {
                retention_time += std::stod(matches[1]) * 60 * 60;
            }
            scan.retention_time = retention_time;

            // Check if we are on the desired region as defined by
            // Grid::Bounds.
            if (retention_time < min_rt) {
                continue;
            }
            // Assuming linearity of the retention time on the mzXML file.
            // We stop searching for the next scan, since we are out of bounds.
            if (retention_time > max_rt) {
                break;
            }
            if (retention_time < raw_data.min_rt) {
                raw_data.min_rt = retention_time;
            }
            if (retention_time > raw_data.max_rt) {
                raw_data.max_rt = retention_time;
            }

            // Fetch the peaks tag.
            tag = XmlReader::read_tag(stream);
            while (tag) {
                if (tag.value().name == "scan" && tag.value().closed) {
                    break;
                }
                if (tag.value().name == "peaks") {
                    auto peak_attributes = tag.value().attributes;

                    // Extract the precision from the peaks tag.
                    if (peak_attributes.find("precision") ==
                        peak_attributes.end()) {
                        return std::nullopt;
                    }
                    int precision = std::stoi(peak_attributes["precision"]);

                    // Extract the byteOrder from the peaks tag. This determines
                    // the endianness in which the data was stored. `network` ==
                    // `big_endian`.
                    if (peak_attributes.find("byteOrder") ==
                        peak_attributes.end()) {
                        return std::nullopt;
                    }
                    auto byte_order = peak_attributes["byteOrder"];
                    auto little_endian = byte_order != "network";

                    // Extract the contentType/pairOrder from the peaks tag and
                    // exit if it is not `m/z-int`. In older versions of
                    // ProteoWizard, the conversion was not validated and the
                    // tag was incorrect. Here we are supporting both versions
                    // for compatibility but we are not trying to be exhaustive.
                    auto content_type_found =
                        peak_attributes.find("contentType") !=
                        peak_attributes.end();
                    auto pair_order_found = peak_attributes.find("pairOrder") !=
                                            peak_attributes.end();
                    if (!content_type_found && !pair_order_found) {
                        return std::nullopt;
                    }
                    std::string pair_order;
                    if (content_type_found) {
                        pair_order = peak_attributes["contentType"];
                    } else {
                        pair_order = peak_attributes["pairOrder"];
                    }

                    // Extract the peaks from the data.
                    auto data = XmlReader::read_data(stream);
                    if (!data) {
                        return std::nullopt;
                    }

                    // Decode the points from the base 64 string.
                    // Initialize Base64 decoder.
                    Base64 decoder(
                        reinterpret_cast<unsigned char *>(&data.value()[0]),
                        precision, little_endian);
                    double intensity_sum = 0;
                    double max_intensity = 0;
                    for (size_t i = 0; i < num_points; ++i) {
                        auto mz = decoder.get_double();
                        auto intensity = decoder.get_double();

                        // We don't need to extract the peaks when we are not
                        // inside the mz bounds or contain no value.
                        if (mz < min_mz || mz > max_mz || intensity == 0) {
                            continue;
                        }
                        if (mz < raw_data.min_mz) {
                            raw_data.min_mz = mz;
                        }
                        if (mz > raw_data.max_mz) {
                            raw_data.max_mz = mz;
                        }

                        if (intensity > max_intensity) {
                            max_intensity = intensity;
                        }
                        intensity_sum += intensity;

                        // NOTE(alex): Not the most efficient way. It would be
                        // better to preallocate but we don't know at this point
                        // how many peaks from peak_count are inside our bounds.
                        scan.mz.push_back(mz);
                        scan.intensity.push_back(intensity);
                    }
                    scan.num_points = scan.mz.size();
                    if (scan.num_points != 0) {
                        raw_data.scans.push_back(scan);
                        raw_data.retention_times.push_back(retention_time);
                        scan.max_intensity = max_intensity;
                        scan.total_intensity = intensity_sum;
                    }
                }
                if (tag.value().name == "precursorMz") {
                    auto precursor_attributes = tag.value().attributes;
                    if (precursor_attributes.find("precursorIntensity") !=
                        precursor_attributes.end()) {
                        scan.precursor_information.intensity = std::stod(
                            precursor_attributes["precursorIntensity"]);
                    }

                    if (precursor_attributes.find("windowWideness") !=
                        precursor_attributes.end()) {
                        scan.precursor_information.window_wideness =
                            std::stod(precursor_attributes["windowWideness"]);
                    }

                    if (precursor_attributes.find("windowWideness") !=
                        precursor_attributes.end()) {
                        scan.precursor_information.window_wideness =
                            std::stod(precursor_attributes["windowWideness"]);
                    }

                    if (precursor_attributes.find("precursorCharge") !=
                        precursor_attributes.end()) {
                        scan.precursor_information.charge =
                            std::stoi(precursor_attributes["precursorCharge"]);
                    }

                    if (precursor_attributes.find("activationMethod") !=
                        precursor_attributes.end()) {
                        if (precursor_attributes["activationMethod"] == "CID") {
                            scan.precursor_information.activation_method =
                                ActivationMethod::CID;
                        } else if (precursor_attributes["activationMethod"] ==
                                   "HCD") {
                            scan.precursor_information.activation_method =
                                ActivationMethod::HCD;
                        } else {
                            scan.precursor_information.activation_method =
                                ActivationMethod::UNKNOWN;
                        }
                    } else {
                        scan.precursor_information.activation_method =
                            ActivationMethod::UNKNOWN;
                    }

                    auto data = XmlReader::read_data(stream);
                    if (!data) {
                        return std::nullopt;
                    }
                    scan.precursor_information.mz = std::stod(data.value());

                    if (precursor_attributes.find("precursorScanNum") !=
                        precursor_attributes.end()) {
                        scan.precursor_information.scan_number =
                            std::stoi(precursor_attributes["precursorScanNum"]);
                    } else {
                        scan.precursor_information.scan_number =
                            previous_ms1_scan_number;
                    }
                }
                tag = XmlReader::read_tag(stream);
            }
        }
    }
    return raw_data;
}

std::optional<std::string> XmlReader::read_data(std::istream &stream) {
    std::string data;
    std::getline(stream, data, '<');
    if (data.empty() || !stream.good()) {
        return std::nullopt;
    }

    // Trim potential whitespace at the beginning of the data string.
    for (size_t i = 0; i < data.size(); ++i) {
        if (!std::isspace(data[i])) {
            if (i != 0) {
                data = data.substr(i);
            }
            break;
        }
    }

    return data;
}

std::optional<XmlReader::Tag> XmlReader::read_tag(std::istream &stream) {
    bool is_closed = false;
    bool reading_content = false;
    auto tag = std::optional<Tag>(Tag{});

    // Store the tag contents in a buffer for further processing.
    std::string buffer;
    while (stream.good()) {
        char c = stream.get();
        if (c == '<') {
            if (stream.peek() == '/') {
                tag->closed = true;
                stream.get();
            } else if (stream.peek() == ' ') {
                return std::nullopt;
            }
            reading_content = true;
        } else if (c == '>') {
            is_closed = true;
            break;
        } else if (reading_content) {
            buffer += c;
        }
    }

    if (!is_closed) {
        return std::nullopt;
    }

    // Tokenize tag contents.
    std::stringstream ss(buffer);
    ss >> tag->name;

    // Find all attributes of this tag and store them on the attribute map.
    std::regex attribute_regex("(\\S+)=\"([^\"]+)\"");
    std::smatch matches;
    while (std::regex_search(buffer, matches, attribute_regex)) {
        if (matches.size() != 3 || matches[1] == "" || matches[2] == "") {
            return std::nullopt;
        }
        tag->attributes[matches[1]] = matches[2];
        buffer = matches.suffix().str();
    }

    return tag;
};

IdentData::IdentData XmlReader::read_mzidentml(std::istream &stream) {
    std::vector<IdentData::SpectrumId> spectrum_ids;
    std::vector<IdentData::Peptide> peptides;
    std::vector<IdentData::DBSequence> db_sequences;
    std::vector<IdentData::ProteinHypothesis> protein_hypotheses;

    // Find all Peptides and DBSequences in the file (SequenceCollection).
    while (stream.good()) {
        auto tag = XmlReader::read_tag(stream);
        if (!tag) {
            continue;
        }
        if (tag.value().name == "SequenceCollection" && tag.value().closed) {
            break;
        }
        if (tag.value().name == "Peptide") {
            auto peptide = IdentData::Peptide{};
            auto attributes = tag.value().attributes;
            peptide.id = attributes["id"];
            while (stream.good()) {
                auto tag = XmlReader::read_tag(stream);
                if (tag.value().name == "Peptide" && tag.value().closed) {
                    peptides.push_back(peptide);
                    break;
                }
                if (tag.value().name == "PeptideSequence" &&
                    !tag.value().closed) {
                    auto data = XmlReader::read_data(stream);
                    if (!data) {
                        break;
                        // FIXME: Throw exception? Return nullopt?
                        // return std::nullopt;
                    }
                    peptide.sequence = data.value();
                }
                if (tag.value().name == "Modification" && !tag.value().closed) {
                    // Save modification info.
                    auto attributes = tag.value().attributes;
                    auto modification = IdentData::PeptideModification{};
                    if (attributes.find("monoisotopicMassDelta") !=
                        attributes.end()) {
                        modification.monoisotopic_mass_delta =
                            std::stod(attributes["monoisotopicMassDelta"]);
                    }
                    if (attributes.find("avgMassDelta") != attributes.end()) {
                        modification.average_mass_delta =
                            std::stod(attributes["avgMassDelta"]);
                    }
                    if (attributes.find("residues") != attributes.end()) {
                        modification.residues = attributes["residues"];
                    }
                    if (attributes.find("location") != attributes.end()) {
                        modification.location =
                            std::stoi(attributes["location"]);
                    }
                    // Find CVParams for this modification..
                    while (stream.good()) {
                        auto tag = XmlReader::read_tag(stream);
                        if (tag.value().name == "cvParam") {
                            auto cv_param = IdentData::CVParam{};
                            auto attributes = tag.value().attributes;
                            cv_param.name = attributes["name"];
                            cv_param.accession = attributes["accession"];
                            cv_param.cv_ref = attributes["cvRef"];
                            if (attributes.find("value") != attributes.end()) {
                                cv_param.value = attributes["value"];
                            }
                            modification.cv_params.push_back(cv_param);
                        }
                        if (tag.value().name == "Modification" &&
                            tag.value().closed) {
                            peptide.modifications.push_back(modification);
                            break;
                        }
                    }
                }
            }
        }
        if (tag.value().name == "DBSequence") {
            auto db_sequence = IdentData::DBSequence{};
            auto attributes = tag.value().attributes;
            db_sequence.id = attributes["id"];
            while (stream.good()) {
                auto tag = XmlReader::read_tag(stream);
                if (tag.value().name == "DBSequence" && tag.value().closed) {
                    db_sequences.push_back(db_sequence);
                    break;
                }
                auto attributes = tag.value().attributes;
                if (tag.value().name == "cvParam" &&
                    attributes["accession"] == "MS:1001088") {
                    db_sequence.value = attributes["value"];
                }
            }
        }
    }
    // Find the PSMs for this data (SpectrumIdentificationResult).
    while (stream.good()) {
        auto tag = XmlReader::read_tag(stream);
        if (!tag) {
            continue;
        }
        if (tag.value().name == "SpectrumIdentificationList" &&
            tag.value().closed) {
            break;
        }
        if (tag.value().name != "SpectrumIdentificationResult") {
            continue;
        }
        auto spectrum_id = IdentData::SpectrumId{};
        bool identification_item_found = false;
        while (stream.good()) {
            tag = XmlReader::read_tag(stream);
            auto attributes = tag.value().attributes;

            // Retention time.
            if (tag.value().name == "cvParam" &&
                attributes["accession"] == "MS:1000894") {
                spectrum_id.retention_time = std::stod(attributes["value"]);
            }

            // Identification item.
            if (tag.value().name == "SpectrumIdentificationItem" &&
                !tag.value().closed) {
                if (identification_item_found &&
                    std::stoi(attributes["rank"]) < spectrum_id.rank) {
                    continue;
                }
                spectrum_id.id = attributes["id"];
                spectrum_id.rank = std::stoi(attributes["rank"]);
                spectrum_id.pass_threshold =
                    attributes["passThreshold"] == "true";
                spectrum_id.peptide_id = attributes["peptide_ref"];
                spectrum_id.charge_state = std::stoi(attributes["chargeState"]);
                spectrum_id.theoretical_mz =
                    std::stod(attributes["calculatedMassToCharge"]);
                spectrum_id.experimental_mz =
                    std::stod(attributes["experimentalMassToCharge"]);
                identification_item_found = true;
            }

            if (tag.value().name == "SpectrumIdentificationResult" &&
                tag.value().closed) {
                spectrum_ids.push_back(spectrum_id);
                break;
            }
        }
    }
    // Find the protein groups for this data (ProteinDetectionList).
    while (stream.good()) {
        auto tag = XmlReader::read_tag(stream);
        if (!tag) {
            continue;
        }
        if (tag.value().name == "ProteinDetectionList" && tag.value().closed) {
            break;
        }
        if (tag.value().name == "ProteinDetectionHypothesis" &&
            !tag.value().closed) {
            auto protein_hypothesis = IdentData::ProteinHypothesis{};
            auto attributes = tag.value().attributes;
            protein_hypothesis.db_sequence_id = attributes["dBSequence_ref"];
            protein_hypothesis.pass_threshold =
                attributes["passThreshold"] == "true";
            while (stream.good()) {
                tag = XmlReader::read_tag(stream);
                if (tag.value().name == "ProteinDetectionHypothesis" &&
                    tag.value().closed) {
                    protein_hypotheses.push_back(protein_hypothesis);
                    break;
                }
                if (tag.value().name == "SpectrumIdentificationItemRef") {
                    auto attributes = tag.value().attributes;
                    protein_hypothesis.spectrum_ids.push_back(
                        attributes["spectrumIdentificationItem_ref"]);
                }
            }
        }
    }
    // Cross link peptide_id per SpectrumId to obtain the original sequence.
    for (auto &ident : spectrum_ids) {
        for (const auto &peptide : peptides) {
            if (peptide.id == ident.peptide_id) {
                ident.sequence = peptide.sequence;
                if (!peptide.modifications.empty()) {
                    ident.modifications = true;
                }
                break;
            }
        }
    }
    return {db_sequences, peptides, spectrum_ids, protein_hypotheses};
}
