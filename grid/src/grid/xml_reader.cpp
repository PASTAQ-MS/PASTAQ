#include <regex>
#include <sstream>

#include "utils/base64.hpp"
#include "xml_reader.hpp"

namespace RawData {
enum Polarity : unsigned char { POSITIVE, NEGATIVE, BOTH };
enum ActivationMethod : unsigned char { CID, HCD };

struct PrecursorInformation {
    // Index for the scan that caused the MSn event.
    uint64_t scan_number;
    // Detected charge for the precursor.
    uint8_t charge;
    // Mass to charge of the MSn event.
    double mz;
    // Intensity of the precursor event.
    double intensity;
    // The activation method for the fragmentation of the MSn event.
    ActivationMethod activation_method;
    // The total isolation window selected for fragmentation in m/z units.
    double window_wideness;
};

struct Scan {
    // Index for this scan.
    uint64_t scan_number;
    // Type of ms_level of this scan (i.e. MS1/MS2/MSn).
    uint64_t ms_level;
    // How many mz-intensity pairs are containd in this scan.
    uint64_t num_points;
    // Retention time in seconds of this scan;
    double retention_time;
    // mz-intensity vectors should have the same size (num_points).
    std::vector<double> mz;
    std::vector<double> intensity;
    // The polarity of the ionization for this scan.
    Polarity polarity;
    // In case this is a MSn scan, the precursor information will be stored
    // here.
    PrecursorInformation *precursor_information;
};

struct RawData {
    // The instrument type.
    Instrument::Type instrument_type;
    // Min/max mass to charge range (m/z).
    double min_mz;
    double max_mz;
    // Min/max retention time range (seconds).
    double min_rt;
    double max_rt;
    // Resolution of MS1/MSn at the reference m/z. In this case the resolution
    // is defined as:
    //
    //     R = reference_mz/fwhm_at_reference_mz
    //
    double resolution_ms1;
    double resolution_msn;
    double reference_mz;

    // TODO(alex): Add more metadata information here.

    // Extracted scans.
    std::vector<::RawData::Scan> scans;
};
}  // namespace RawData

// TODO: Test this function and ensure the results are the same as
// using the other reader.
std::optional<RawData::RawData> read_mzxml(
    std::istream &stream, double min_mz, double max_mz, double min_rt,
    double max_rt, Instrument::Type instrument_type, double resolution_ms1,
    double resolution_msn, double reference_mz, RawData::Polarity polarity) {
    RawData::RawData raw_data{instrument_type,
                              std::numeric_limits<double>::infinity(),
                              -std::numeric_limits<double>::infinity(),
                              std::numeric_limits<double>::infinity(),
                              -std::numeric_limits<double>::infinity(),
                              resolution_ms1,
                              resolution_msn,
                              reference_mz,
                              {}};
    // TODO(alex): Can we automatically detect the instrument type and set
    // resolution from the header?
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
            if (scan_attributes.find("polarity") == scan_attributes.end()) {
                return std::nullopt;
            }
            if (scan_attributes["polarity"] == "+") {
                scan.polarity = RawData::Polarity::POSITIVE;
            } else if (scan_attributes["polarity"] == "-") {
                scan.polarity = RawData::Polarity::NEGATIVE;
            }
            if (scan.polarity != RawData::Polarity::BOTH &&
                scan.polarity != polarity) {
                return std::nullopt;
            }

            // Find MS level.
            if (scan_attributes.find("msLevel") == scan_attributes.end()) {
                return std::nullopt;
            }
            // NOTE(alex): We are only interested on ms1 scans for now.
            if (scan_attributes["msLevel"] != "1") {
                continue;
            }
            size_t ms_level = std::stoi(scan_attributes["msLevel"]);
            scan.ms_level = ms_level;

            // Find the number of m/z-intensity pairs in the scan.
            if (scan_attributes.find("peaksCount") == scan_attributes.end()) {
                return std::nullopt;
            }
            size_t num_points = std::stoi(scan_attributes["peaksCount"]);

            // Extract the retention time.
            if (scan_attributes.find("retentionTime") ==
                scan_attributes.end()) {
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
            // We stop searching for the scan, since we are out of bounds.
            if (retention_time > max_rt) {
                return std::nullopt;
            }

            // TODO(alex): Extract the PrecursorInformation in case it is a MSn
            // scan.

            // Fetch the peaks tag.
            tag = XmlReader::read_tag(stream);
            while (tag) {
                if (tag.value().name == "scan" && tag.value().closed) {
                    return std::nullopt;
                }
                if (tag.value().name == "peaks") {
                    break;
                }
                tag = XmlReader::read_tag(stream);
            }
            auto peak_attributes = tag.value().attributes;

            // Extract the precision from the peaks tag.
            if (peak_attributes.find("precision") == peak_attributes.end()) {
                return std::nullopt;
            }
            int precision = std::stoi(peak_attributes["precision"]);

            // Extract the byteOrder from the peaks tag. This determines the
            // endianness in which the data was stored. `network` ==
            // `big_endian`.
            if (peak_attributes.find("byteOrder") == peak_attributes.end()) {
                return std::nullopt;
            }
            auto byte_order = peak_attributes["byteOrder"];
            auto little_endian = byte_order != "network";

            // Extract the contentType/pairOrder from the peaks tag and exit
            // if it is not `m/z-int`. In older versions of ProteoWizard,
            // the conversion was not validated and the tag was incorrect.
            // Here we are supporting both versions for compatibility but we
            // are not trying to be exhaustive.
            auto content_type_found =
                peak_attributes.find("contentType") != peak_attributes.end();
            auto pair_order_found =
                peak_attributes.find("pairOrder") != peak_attributes.end();
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
            Base64 decoder(reinterpret_cast<unsigned char *>(&data.value()[0]),
                           precision, little_endian);
            for (size_t i = 0; i < num_points; ++i) {
                auto mz = decoder.get_double();
                auto intensity = decoder.get_double();

                // We don't need to extract the peaks when we are not inside
                // the mz bounds or contain no value.
                if (mz < min_mz || mz > max_mz || intensity == 0) {
                    continue;
                }

                // NOTE(alex): Not the most efficient way. It would be better to
                // preallocate but we don't know at this point how many
                // peaks from peak_count are inside our bounds.
                scan.mz.push_back(mz);
                scan.intensity.push_back(intensity);
            }
            // TODO(alex): ASSERT scan.mz.size() == scan.intensity.size()
            scan.num_points = scan.mz.size();
            raw_data.scans.push_back(scan);
        }
    }
    return raw_data;
}

// Read the next mz1 scan from the stream.
std::optional<std::vector<Grid::Point>> XmlReader::read_next_scan(
    std::istream &stream, const Grid::Parameters &parameters) {
    // TODO(alex): This function does not need Grid::Parameters, only the
    // Grid::Bounds.
    while (stream.good()) {
        auto tag = XmlReader::read_tag(stream);
        if (!tag) {
            continue;
        }
        if (tag.value().name == "scan" && !tag.value().closed) {
            auto scan_attributes = tag.value().attributes;
            // We are only interested on ms1 scans for now.
            if (scan_attributes["msLevel"] != "1") {
                continue;
            }
            // Extract the peak count.
            if (scan_attributes.find("peaksCount") == scan_attributes.end()) {
                return std::nullopt;
            }
            size_t peaks_count = std::stoi(scan_attributes["peaksCount"]);

            // Extract the retention time.
            if (scan_attributes.find("retentionTime") ==
                scan_attributes.end()) {
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

            // Check if we are on the desired region as defined by
            // Grid::Bounds.
            if (retention_time < parameters.bounds.min_rt) {
                continue;
            }
            // Assuming linearity of the retention time on the mzXML file.
            // We stop searching for the scan, since we are out of bounds.
            if (retention_time > parameters.bounds.max_rt) {
                return std::nullopt;
            }

            // Fetch the peaks tag.
            tag = XmlReader::read_tag(stream);
            while (tag) {
                if (tag.value().name == "scan" && tag.value().closed) {
                    return std::nullopt;
                }
                if (tag.value().name == "peaks") {
                    break;
                }
                tag = XmlReader::read_tag(stream);
            }
            auto peak_attributes = tag.value().attributes;

            // Extract the precision from the peaks tag.
            if (peak_attributes.find("precision") == peak_attributes.end()) {
                return std::nullopt;
            }
            int precision = std::stoi(peak_attributes["precision"]);

            // Extract the byteOrder from the peaks tag. This determines the
            // endianness in which the data was stored. `network` ==
            // `big_endian`.
            if (peak_attributes.find("byteOrder") == peak_attributes.end()) {
                return std::nullopt;
            }
            auto byte_order = peak_attributes["byteOrder"];
            auto little_endian = byte_order != "network";

            // Extract the contentType/pairOrder from the peaks tag and exit
            // if it is not `m/z-int`. In older versions of ProteoWizard,
            // the conversion was not validated and the tag was incorrect.
            // Here we are supporting both versions for compatibility but we
            // are not trying to be exhaustive.
            auto content_type_found =
                peak_attributes.find("contentType") != peak_attributes.end();
            auto pair_order_found =
                peak_attributes.find("pairOrder") != peak_attributes.end();
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
            std::vector<Grid::Point> points;

            // Initialize Base64 decoder.
            Base64 decoder(reinterpret_cast<unsigned char *>(&data.value()[0]),
                           precision, little_endian);
            for (size_t i = 0; i < peaks_count; ++i) {
                auto mz = decoder.get_double();
                auto intensity = decoder.get_double();

                // We don't need to extract the peaks when we are not inside
                // the mz bounds or contain no value.
                if (mz < parameters.bounds.min_mz ||
                    mz > parameters.bounds.max_mz || intensity == 0) {
                    continue;
                }

                // Not the most efficient way. It would be better to
                // preallocate but we don't know at this point how many
                // peaks from peak_count are inside our bounds.
                points.push_back({mz, retention_time, intensity});
            }

            return points;
        }
    }
    return std::nullopt;
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
