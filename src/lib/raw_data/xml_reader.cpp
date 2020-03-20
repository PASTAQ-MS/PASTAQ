#include <regex>
#include <sstream>

#include "utils/base64.hpp"
#include "xml_reader.hpp"

RawData::Scan parse_mzxml_scan(std::istream &stream,
                               std::optional<XmlReader::Tag> &tag,
                               double min_mz, double max_mz, double min_rt,
                               double max_rt, Polarity::Type polarity,
                               size_t ms_level) {
    RawData::Scan scan;
    uint64_t precursor_id = 0;
    scan.precursor_information.scan_number = 0;
    auto scan_attributes = tag.value().attributes;

    // Find scan number.
    if (scan_attributes.find("num") == scan_attributes.end()) {
        return {};
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
            return {};
        }
    }

    // Find MS level.
    if (scan_attributes.find("msLevel") == scan_attributes.end()) {
        return {};
    }
    size_t scan_ms_level = std::stoi(scan_attributes["msLevel"]);
    scan.ms_level = scan_ms_level;

    // Fill up the rest of the scan information.
    if (scan_ms_level == ms_level) {
        // Find the number of m/z-intensity pairs in the scan.
        if (scan_attributes.find("peaksCount") == scan_attributes.end()) {
            return {};
        }
        size_t num_points = std::stoi(scan_attributes["peaksCount"]);

        // Extract the retention time.
        if (scan_attributes.find("retentionTime") == scan_attributes.end()) {
            // NOTE(alex): On the spec, the retention time attribute is
            // optional, however, we do require it.
            return {};
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
            return {};
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
        if (retention_time < min_rt) {
            return {};
        }
        // Assuming linearity of the retention time on the mzXML file.
        // TODO: We should stop searching for the next scan, since we are out of
        // bounds.
        if (retention_time > max_rt) {
            return {};
        }

        // Fetch the next tag. We are interested in the contents of this scan
        // tag: precursorMz and peaks.
        auto next_tag = XmlReader::read_tag(stream);
        while (next_tag) {
            if (next_tag.value().name == "scan" && next_tag.value().closed) {
                break;
            }
            // We are not interested in subscans here. Continue until we find
            // the corresponding output.
            // FIXME: This only works with ONE subscan! ms1->ms2, if we are in
            // ms1 and our ms2 contain subscans, it will FAIL. We need to
            // recursively check child scans.
            if (next_tag.value().name == "scan" && !next_tag.value().closed) {
                next_tag = XmlReader::read_tag(stream);
                while (next_tag) {
                    if (next_tag.value().name == "scan" &&
                        next_tag.value().closed) {
                        break;
                    }
                    next_tag = XmlReader::read_tag(stream);
                }
            }
            if (next_tag.value().name == "peaks") {
                auto peak_attributes = next_tag.value().attributes;

                // Extract the precision from the peaks tag.
                if (peak_attributes.find("precision") ==
                    peak_attributes.end()) {
                    return {};
                }
                int precision = std::stoi(peak_attributes["precision"]);

                // Extract the byteOrder from the peaks tag. This determines
                // the endianness in which the data was stored. `network` ==
                // `big_endian`.
                if (peak_attributes.find("byteOrder") ==
                    peak_attributes.end()) {
                    return {};
                }
                auto byte_order = peak_attributes["byteOrder"];
                auto little_endian = byte_order != "network";

                // Extract the contentType/pairOrder from the peaks tag and
                // exit if it is not `m/z-int`. In older versions of
                // ProteoWizard, the conversion was not validated and the
                // tag was incorrect. Here we are supporting both versions
                // for compatibility but we are not trying to be exhaustive.
                auto content_type_found = peak_attributes.find("contentType") !=
                                          peak_attributes.end();
                auto pair_order_found =
                    peak_attributes.find("pairOrder") != peak_attributes.end();
                if (!content_type_found && !pair_order_found) {
                    return {};
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
                    return {};
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
                scan.max_intensity = max_intensity;
                scan.total_intensity = intensity_sum;
            }
            if (next_tag.value().name == "precursorMz") {
                auto precursor_attributes = next_tag.value().attributes;
                if (precursor_attributes.find("precursorIntensity") !=
                    precursor_attributes.end()) {
                    scan.precursor_information.intensity =
                        std::stod(precursor_attributes["precursorIntensity"]);
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
                    return {};
                }
                scan.precursor_information.mz = std::stod(data.value());

                if (precursor_attributes.find("precursorScanNum") !=
                    precursor_attributes.end()) {
                    scan.precursor_information.scan_number =
                        std::stoi(precursor_attributes["precursorScanNum"]);
                }
            }
            next_tag = XmlReader::read_tag(stream);
        }
        return scan;
    }

    // If we are interested in a scan with a mz_level != 1, we need to find the
    // precursor to which it belongs.
    if (scan_ms_level == ms_level - 1) {
        precursor_id = scan.scan_number;
        auto next_tag = XmlReader::read_tag(stream);
        while (next_tag) {
            if (next_tag.value().name == "scan" && !next_tag.value().closed) {
                auto child_scan =
                    parse_mzxml_scan(stream, next_tag, min_mz, max_mz, min_rt,
                                     max_rt, polarity, ms_level);
                child_scan.precursor_information.scan_number = precursor_id;
                if (child_scan.precursor_information.scan_number !=
                        precursor_id &&
                    child_scan.precursor_information.scan_number == 0) {
                    child_scan.precursor_information.scan_number = precursor_id;
                }
                return child_scan;
            }
            next_tag = XmlReader::read_tag(stream);
        }
    }

    return {};
}

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
    while (stream.good()) {
        auto tag = XmlReader::read_tag(stream);
        if (!tag) {
            continue;
        }
        if (tag.value().name == "scan" && !tag.value().closed) {
            auto scan = parse_mzxml_scan(stream, tag, min_mz, max_mz, min_rt,
                                         max_rt, polarity, ms_level);
            if (scan.num_points != 0) {
                raw_data.scans.push_back(scan);
                raw_data.retention_times.push_back(scan.retention_time);
                if (scan.retention_time < raw_data.min_rt) {
                    raw_data.min_rt = scan.retention_time;
                }
                if (scan.retention_time > raw_data.max_rt) {
                    raw_data.max_rt = scan.retention_time;
                }
                if (scan.mz[0] < raw_data.min_mz) {
                    raw_data.min_mz = scan.mz[0];
                }
                if (scan.mz[scan.mz.size() - 1] > raw_data.max_mz) {
                    raw_data.max_mz = scan.mz[scan.mz.size() - 1];
                }
            }
        }
    }
    return raw_data;
}

std::optional<RawData::RawData> XmlReader::read_mzml(
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
    while (stream.good()) {
        auto tag = XmlReader::read_tag(stream);
        if (!tag) {
            continue;
        }
        if (tag.value().name == "spectrumList" && tag.value().closed) {
            break;
        }
        RawData::Scan scan = {};
        if (tag.value().name == "spectrum" && !tag.value().closed) {
            // Parse the contents and metadata of this spectrum.
            scan.precursor_information.scan_number = 0;
            auto scan_attributes = tag.value().attributes;

            // NOTE: In the mzML spec, the native scan number is described on
            // the "id" attribute, and can contain more information than
            // required for just an integer identifer. Moreover, it looks like,
            // at least for Orbitrap data, the scan numbers are non-zero
            // consecutive integers. For the sake of time, I'm just assuming
            // here that this assumption is the same for all formats, but should
            // probably find a more robust way of doing this.
            scan.scan_number = std::stoi(scan_attributes["index"]) + 1;
            std::vector<bool> filter_points;
            std::vector<double> mzs;
            std::vector<double> intensities;
            while (stream.good()) {
                auto tag = XmlReader::read_tag(stream);
                if (!tag) {
                    break;
                }
                if (tag.value().name == "spectrum" && tag.value().closed) {
                    break;
                }

                if (tag.value().name == "cvParam") {
                    auto cv_attributes = tag.value().attributes;
                    auto accession = cv_attributes["accession"];

                    // This scan is ms_level 1
                    if (accession == "MS:1000579") {
                        scan.ms_level = 1;
                    }

                    // MS level a multi-level MSn experiment.
                    if (accession == "MS:1000511") {
                        size_t scan_ms_level =
                            std::stoi(cv_attributes["value"]);
                        scan.ms_level = scan_ms_level;
                    }

                    // Polarity.
                    if (accession == "MS:1000130") {
                        if (polarity != Polarity::BOTH &&
                            Polarity::POSITIVE != polarity) {
                            continue;
                        }
                        scan.polarity = Polarity::POSITIVE;
                    }
                    if (accession == "MS:1000129") {
                        if (polarity != Polarity::BOTH &&
                            Polarity::NEGATIVE != polarity) {
                            continue;
                        }
                        scan.polarity = Polarity::NEGATIVE;
                    }

                    // Retention time.
                    if (accession == "MS:1000016") {
                        scan.retention_time = std::stod(cv_attributes["value"]);
                        // Retention time is store in seconds. Make sure it is
                        // the right unit. If the unit accession was
                        // "UO:0000010" it would be in seconds, so no action is
                        // required.
                        if (cv_attributes["unitAccession"] == "UO:0000031") {
                            scan.retention_time *= 60.0;
                        }
                        if (scan.retention_time < min_rt ||
                            scan.retention_time > max_rt) {
                            scan = {};
                            break;
                        }
                    }
                }

                if (tag.value().name == "precursor") {
                    // TODO: Get precursor id from "spectrumRef" attribute.
                    scan.precursor_information.scan_number = 0;
                    scan.precursor_information.charge = 0;
                    scan.precursor_information.mz = 0.0;
                    scan.precursor_information.window_wideness = 0.0;
                    scan.precursor_information.intensity = 0.0;
                    scan.precursor_information.activation_method =
                        ActivationMethod::UNKNOWN;
                    while (stream.good()) {
                        auto tag = XmlReader::read_tag(stream);
                        if (!tag) {
                            break;
                        }
                        if (tag.value().name == "precursor" &&
                            tag.value().closed) {
                            break;
                        }
                        if (tag.value().name == "cvParam") {
                            auto cv_attributes = tag.value().attributes;
                            auto accession = cv_attributes["accession"];
                            // Isolation window.
                            if (accession == "MS:1000827") {
                                scan.precursor_information.mz =
                                    std::stod(cv_attributes["value"]);
                            }
                            if (accession == "MS:1000828") {
                                scan.precursor_information.window_wideness +=
                                    std::stod(cv_attributes["value"]);
                            }
                            if (accession == "MS:1000829") {
                                scan.precursor_information.window_wideness +=
                                    std::stod(cv_attributes["value"]);
                            }
                            // Charge state.
                            if (accession == "MS:1000041") {
                                scan.precursor_information.charge =
                                    std::stoi(cv_attributes["value"]);
                            }
                            if (accession == "MS:1000042") {
                                scan.precursor_information.intensity =
                                    std::stod(cv_attributes["value"]);
                            }
                            // Activation method.
                            if (accession == "MS:1000422") {
                                scan.precursor_information.activation_method =
                                    ActivationMethod::HCD;
                            }
                        }
                    }
                }

                if (tag.value().name == "binaryDataArray") {
                    // precision can be: 64 or 32 (bits).
                    int precision = 0;
                    // mz: 0, intensity: 1
                    int type = -1;
                    std::optional<std::string> data;
                    size_t num_points =
                        std::stoi(tag.value().attributes["encodedLength"]);
                    while (stream.good()) {
                        auto tag = XmlReader::read_tag(stream);
                        if (!tag) {
                            break;
                        }
                        if (tag.value().name == "binaryDataArray" &&
                            tag.value().closed) {
                            break;
                        }
                        if (tag.value().name == "cvParam") {
                            auto cv_attributes = tag.value().attributes;
                            auto accession = cv_attributes["accession"];
                            // Precision.
                            if (accession == "MS:1000523") {
                                precision = 64;
                            }
                            if (accession == "MS:1000521") {
                                precision = 32;
                            }
                            // Type of vector.
                            if (accession == "MS:1000514") {
                                type = 0;
                            }
                            if (accession == "MS:1000515") {
                                type = 1;
                            }
                        }
                        if (tag.value().name == "binary" &&
                            !tag.value().closed) {
                            data = XmlReader::read_data(stream);
                        }
                    }
                    if (data) {
                        Base64 decoder(
                            reinterpret_cast<unsigned char *>(&data.value()[0]),
                            precision, true);
                        num_points = num_points / (precision / 8) / 4 * 3;
                        if (type == 0) {  // mz
                            mzs = std::vector<double>(num_points);
                        }
                        if (type == 1) {  // intensity
                            intensities = std::vector<double>(num_points);
                        }
                        if (filter_points.empty()) {
                            filter_points =
                                std::vector<bool>(num_points, false);
                        }
                        for (size_t i = 0; i < num_points; ++i) {
                            auto value = decoder.get_double();
                            if (type == 0) {  // mz
                                mzs[i] = value;
                                if (value < min_mz || value > max_mz) {
                                    filter_points[i] = true;
                                }
                            }
                            if (type == 1) {  // intensity
                                intensities[i] = value;
                                if (value == 0.0) {
                                    filter_points[i] = true;
                                }
                            }
                        }
                    }
                }
            }

            // Filter mzs not in range and intensity == 0 scans and calculate
            // max_intensity and total_intensity.
            double intensity_sum = 0;
            double max_intensity = 0;
            for (size_t i = 0; i < filter_points.size(); ++i) {
                if (filter_points[i]) {
                    continue;
                }
                scan.mz.push_back(mzs[i]);
                scan.intensity.push_back(intensities[i]);
                if (intensities[i] > max_intensity) {
                    max_intensity = intensities[i];
                }
                intensity_sum += intensities[i];
            }
            scan.num_points = scan.mz.size();
            scan.max_intensity = max_intensity;
            scan.total_intensity = intensity_sum;

            // TODO: Assert that mz.size() == intenstiy.size()
            if (scan.ms_level == 0 || scan.ms_level != ms_level) {
                scan = {};
            }
        }

        // Update RawData.
        if (scan.num_points != 0) {
            if (scan.retention_time < min_rt) {
                continue;
            }
            if (scan.retention_time > max_rt) {
                break;
            }
            raw_data.scans.push_back(scan);
            raw_data.retention_times.push_back(scan.retention_time);
            if (scan.retention_time < raw_data.min_rt) {
                raw_data.min_rt = scan.retention_time;
            }
            if (scan.retention_time > raw_data.max_rt) {
                raw_data.max_rt = scan.retention_time;
            }
            if (scan.mz[0] < raw_data.min_mz) {
                raw_data.min_mz = scan.mz[0];
            }
            if (scan.mz[scan.mz.size() - 1] > raw_data.max_mz) {
                raw_data.max_mz = scan.mz[scan.mz.size() - 1];
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
            // Read the previous character to check for self closing tag.
            stream.seekg(-2, stream.cur);
            if (stream.good()) {
                if (stream.get() == '/') {
                    tag->closed = true;
                };
                stream.get();
            }
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

IdentData::IdentData XmlReader::read_mzidentml(std::istream &stream,
                                               bool ignore_decoy,
                                               bool require_threshold,
                                               bool max_rank_only) {
    IdentData::IdentData ident_data = {};

    // Find the DBSequences, Peptides and PeptideEvidence in the
    // SequenceCollection tag.
    while (stream.good()) {
        auto tag = XmlReader::read_tag(stream);
        if (!tag) {
            continue;
        }
        if (tag.value().name == "SequenceCollection" && tag.value().closed) {
            break;
        }
        if (tag.value().name == "DBSequence") {
            auto attributes = tag.value().attributes;
            if (attributes.empty()) {
                continue;
            }
            IdentData::DBSequence db_sequence = {};
            db_sequence.id = attributes["id"];
            db_sequence.accession = attributes["accession"];
            db_sequence.db_reference = attributes["searchDatabase_ref"];
            if (!tag.value().closed) {
                // Check if the DBSequence contains the protein description as a
                // cvParam tag.
                while (stream.good()) {
                    auto tag = XmlReader::read_tag(stream);
                    if (!tag) {
                        continue;
                    }
                    if (tag.value().name == "DBSequence" &&
                        tag.value().closed) {
                        break;
                    }
                    auto attributes = tag.value().attributes;
                    if (tag.value().name == "cvParam" &&
                        attributes["accession"] == "MS:1001088") {
                        db_sequence.description = attributes["value"];
                    }
                }
            }
            ident_data.db_sequences.push_back(db_sequence);
        } else if (tag.value().name == "Peptide") {
            IdentData::Peptide peptide = {};
            auto attributes = tag.value().attributes;
            peptide.id = attributes["id"];
            // Find peptide sequence and modifications.
            while (stream.good()) {
                auto tag = XmlReader::read_tag(stream);
                if (!tag) {
                    continue;
                }
                if (tag.value().name == "Peptide" && tag.value().closed) {
                    break;
                }
                if (tag.value().name == "PeptideSequence" &&
                    !tag.value().closed) {
                    auto data = XmlReader::read_data(stream);
                    if (!data) {
                        return {};
                        break;
                    }
                    peptide.sequence = data.value();
                }
                // Search peptide modifications.
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
                    } else {
                        modification.location = -1;
                    }
                    // Find identification information for this modification.
                    while (stream.good()) {
                        auto tag = XmlReader::read_tag(stream);
                        if (!tag) {
                            continue;
                        }
                        if (tag.value().name == "Modification" &&
                            tag.value().closed) {
                            peptide.modifications.push_back(modification);
                            break;
                        }
                        if (tag.value().name == "cvParam") {
                            auto attributes = tag.value().attributes;
                            modification.id.push_back(attributes["accession"] +
                                                      "|" + attributes["name"]);
                        }
                    }
                    peptide.modifications.push_back(modification);
                } else if (tag.value().name == "SubstitutionModification") {
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
                    } else {
                        modification.location = -1;
                    }
                    modification.id.push_back(
                        "SUBSTITUTION|" + attributes["originalResidue"] + "->" +
                        attributes["replacementResidue"]);
                    peptide.modifications.push_back(modification);
                }
            }
            ident_data.peptides.push_back(peptide);
        } else if (tag.value().name == "PeptideEvidence") {
            IdentData::PeptideEvidence peptide_evidence;
            auto attributes = tag.value().attributes;
            peptide_evidence.id = attributes["id"];
            peptide_evidence.db_sequence_id = attributes["dBSequence_ref"];
            peptide_evidence.peptide_id = attributes["peptide_ref"];
            peptide_evidence.decoy = false;
            if (attributes.find("isDecoy") != attributes.end()) {
                peptide_evidence.decoy = attributes["isDecoy"] == "true";
            }
            if (ignore_decoy && peptide_evidence.decoy) {
                continue;
            }
            ident_data.peptide_evidence.push_back(peptide_evidence);
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
        // Find the next SpectrumIdentificationResult.
        if (tag.value().name != "SpectrumIdentificationResult") {
            continue;
        }

        // Record all SpectrumIdentificationItems for this result.
        std::vector<IdentData::SpectrumMatch> spectrum_matches;
        double retention_time = 0.0;
        while (stream.good()) {
            tag = XmlReader::read_tag(stream);
            if (!tag) {
                continue;
            }
            if (tag.value().name == "SpectrumIdentificationResult" &&
                tag.value().closed) {
                break;
            }
            auto attributes = tag.value().attributes;

            if (tag.value().name == "cvParam") {
                // Retention time or scan start time.
                if (attributes["accession"] == "MS:1000894" ||
                    attributes["accession"] == "MS:1000016") {
                    retention_time = std::stod(attributes["value"]);
                    // If the retention time is in minutes, we convert it back
                    // to seconds.
                    if (attributes["unitAccession"] == "UO:0000031") {
                        retention_time *= 60.0;
                    }
                }
            }

            // Identification item.
            if (tag.value().name == "SpectrumIdentificationItem" &&
                !tag.value().closed) {
                IdentData::SpectrumMatch spectrum_match = {};
                spectrum_match.id = attributes["id"];
                spectrum_match.pass_threshold =
                    attributes["passThreshold"] == "true";
                if (require_threshold && !spectrum_match.pass_threshold) {
                    continue;
                }
                spectrum_match.match_id = attributes["peptide_ref"];
                spectrum_match.charge_state =
                    std::stoi(attributes["chargeState"]);
                spectrum_match.experimental_mz =
                    std::stod(attributes["experimentalMassToCharge"]);
                spectrum_match.retention_time = 0;
                spectrum_match.rank = std::stoi(attributes["rank"]);
                // Might be optional according to the mzIdentML v1.2.0 spec.
                if (attributes.find("calculatedMassToCharge") !=
                    attributes.end()) {
                    spectrum_match.theoretical_mz =
                        std::stod(attributes["calculatedMassToCharge"]);
                } else {
                    spectrum_match.theoretical_mz = 0.0;
                }

                // Try to extract identification scores.
                while (stream.good()) {
                    tag = XmlReader::read_tag(stream);
                    if (!tag) {
                        continue;
                    }
                    if (tag.value().name == "SpectrumIdentificationItem" &&
                        tag.value().closed) {
                        break;
                    }
                    if (tag.value().name == "cvParam") {
                        auto attributes = tag.value().attributes;
                        // Comet.
                        if (attributes["accession"] == "MS:1002252") {
                            spectrum_match.score_comet_xcor =
                                std::stod(attributes["value"]);
                        }
                    }
                }
                spectrum_matches.push_back(spectrum_match);
            }
        }

        if (max_rank_only) {
            IdentData::SpectrumMatch selected_spectrum;
            // Update retention time on the provisional spectrum_matches list
            // and find the maximum rank spectrum. The rank is in descending
            // order of importance, thus rank 1 is the maximum, and bigger
            // numbers are worse.
            for (size_t i = 0; i < spectrum_matches.size(); ++i) {
                auto &spectrum_match = spectrum_matches[i];
                spectrum_match.retention_time = retention_time;
                if (i == 0 || spectrum_match.rank < selected_spectrum.rank) {
                    selected_spectrum = spectrum_match;
                }
            }
            ident_data.spectrum_matches.push_back(selected_spectrum);
        } else {
            // Update retention time on the provisional spectrum_matches list
            // and push each element to the list of PSM.
            for (auto &spectrum_match : spectrum_matches) {
                spectrum_match.retention_time = retention_time;
                ident_data.spectrum_matches.push_back(spectrum_match);
            }
        }
    }
    return ident_data;
}
