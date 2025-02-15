#include <zlib.h>
#include <regex>
#include <sstream>
#include <cctype>
#include <numeric>
#include <chrono> // For timing

#include "utils/base64.hpp"
#include "utils/compression.hpp"
#include "xml_reader.hpp"

#include "mzParser.h"

#include "MSToolkitTypes.h"
#include "MSReader.h"
#include "MSObject.h"
#include "Spectrum.h"

using namespace mzParser;
using namespace std;
using namespace MSToolkit;

ActivationMethod::Type mapActivationMethod(MSActivation activation) {
    switch (activation) {
        case mstCID:   return ActivationMethod::CID;
        case mstECD:   return ActivationMethod::ECD;
        case mstETD:   return ActivationMethod::ETD;
        case mstETDSA: return ActivationMethod::ETDSA;
        case mstPQD:   return ActivationMethod::PQD;
        case mstHCD:   return ActivationMethod::HCD;
        case mstIRMPD: return ActivationMethod::IRMPD;
        case mstSID:   return ActivationMethod::SID;
        case mstNA:    return ActivationMethod::NA;
        default:       return ActivationMethod::UNKNOWN;
    }
}
// Read a scan using mstoolkit and store it in the RawData::Scan structure
//RawData::Scan read_mzxml_scan(const BasicSpectrum& spectrum) {
RawData::Scan read_scan(MSToolkit::Spectrum& spectrum) {
        
    RawData::Scan scan;

    // Set scan number and MS level
    // scan.scan_number = static_cast<uint64_t>(spectrum.getScanNum());
    // scan.ms_level = static_cast<uint64_t>(spectrum.getMSLevel());

    scan.scan_number = spectrum.getScanNumber();
    scan.ms_level = spectrum.getMsLevel();

    // Get retention time
    scan.retention_time = spectrum.getRTime();

    // Get spectrum data points (mz and intensity)
    size_t num_points = spectrum.size();
    scan.num_points = static_cast<uint64_t>(num_points);
    scan.mz.reserve(num_points);
    scan.intensity.reserve(num_points);

    for (size_t j=0;j<spectrum.size();j++) {
        scan.mz.push_back(spectrum.at(j).mz);
        scan.intensity.push_back(spectrum.at(j).intensity);
    }

    // Calculate max and total intensity
    scan.max_intensity = *std::max_element(scan.intensity.begin(), scan.intensity.end());
    scan.total_intensity = std::accumulate(scan.intensity.begin(), scan.intensity.end(), 0.0);

    // Determine polarity
    // Not sure how to do this from the spectrum data structure 
    // https://github.com/mhoopmann/mstoolkit/blob/master/include/Spectrum.h
    // BasicSpectrum: scan.polarity = spectrum.getPositiveScan() ? Polarity::POSITIVE : Polarity::NEGATIVE;

    MSToolkit::MSPrecursorInfo mspinfo = spectrum.getPrecursor();

    // Get precursor information (for MSn scans)
    if (scan.ms_level > 1) {
        scan.precursor_information.mz = mspinfo.mz;
        // It seems that spectrum.getPrecursorIntensity() is not implemented in BasicSpectrum
        // scan.precursor_information.intensity = spectrum.getPrecursorIntensity();
        scan.precursor_information.charge = mspinfo.charge;
        scan.precursor_information.scan_number = mspinfo.precursorScanNumber;

        // The activation method for the fragmentation of the MSn event.
        // ActivationMethod::Type activation_method;
        scan.precursor_information.activation_method = mapActivationMethod(mspinfo.activation);
        
        // The total isolation window selected for fragmentation in m/z units.
        // double window_wideness = isoOffsetUpper + isoOffsetLower;

          //       can we extimate window_wideness from the following structure?
          //           struct MSPrecursorInfo {
          // double mz=0;
          // double monoMz=0;
          // double isoMz=0;
          // int charge=0;
          // MSActivation activation=mstNA;
          // int precursorScanNumber=0;
          // double isoOffsetLower=0;
          // double isoOffsetUpper=0;
          //   };
    }

    return scan;

}

// Read a scan using mstoolkit and store it in the RawData::Scan structure
//RawData::Scan read_mzxml_scan(const BasicSpectrum& spectrum) {
RawData::Scan read_mzxml_scan(BasicSpectrum& spectrum) {
    //mstoolkit functions
    
    // BasicSpectrum spe;
    // MzParser sax(&s)
    
    RawData::Scan scan;

    // Set scan number and MS level
    // scan.scan_number = static_cast<uint64_t>(spectrum.getScanNum());
    // scan.ms_level = static_cast<uint64_t>(spectrum.getMSLevel());

    scan.scan_number = spectrum.getScanNum();
    scan.ms_level = spectrum.getMSLevel();

    // Get retention time in seconds (convert if stored in minutes)
    scan.retention_time = spectrum.getRTime() * 60.0;

    // Get spectrum data points (mz and intensity)
    size_t num_points = spectrum.size();
    scan.num_points = static_cast<uint64_t>(num_points);
    scan.mz.reserve(num_points);
    scan.intensity.reserve(num_points);

    for (size_t i = 0; i < num_points; ++i) {
        specDP dp = spectrum[i]; // Access data point
        scan.mz.push_back(dp.mz);
        scan.intensity.push_back(dp.intensity);
    }

    // Calculate max and total intensity
    scan.max_intensity = *std::max_element(scan.intensity.begin(), scan.intensity.end());
    scan.total_intensity = std::accumulate(scan.intensity.begin(), scan.intensity.end(), 0.0);

    // Determine polarity
    // from raw_data.hpp: "This parameter is meant to be used as a filter of the scans that are going to be
    // read from the raw data." Probably it is not needed anymore as raw data is read using mstoolkit
    scan.polarity = spectrum.getPositiveScan() ? Polarity::POSITIVE : Polarity::NEGATIVE;

    // Get precursor information (for MSn scans)
    if (scan.ms_level > 1 && spectrum.getPrecursorIonCount() > 0) {
        scan.precursor_information.mz = spectrum.getPrecursorMZ();
        // It seems that spectrum.getPrecursorIntensity() is not implemented in BasicSpectrum
        // scan.precursor_information.intensity = spectrum.getPrecursorIntensity();
        scan.precursor_information.charge = spectrum.getPrecursorCharge();
        scan.precursor_information.scan_number = spectrum.getPrecursorScanNum();

        // missing
            // The activation method for the fragmentation of the MSn event.
        // ActivationMethod::Type activation_method;
            // The total isolation window selected for fragmentation in m/z units.
        // double window_wideness;
    }

    return scan;

}

// Read an entire mzxml file into the RawData::RawData data structure filtering usin mstoolkit libraries1
// std::optional<RawData::RawMSData> XmlReader::read_msgdata(
//     std::string &input_file, double min_mz, double max_mz, double min_rt,
//     double max_rt, Instrument::Type instrument_type, double resolution_ms1,
//     double resolution_msn, double reference_mz, Polarity::Type polarity,
//     size_t ms_level) {

//     BasicSpectrum s;
//     // mz file handler
//     MzParser mzfh(&s);
    
//     RawData::RawMSData raw_data = {};
//     raw_data.instrument_type = instrument_type;
//     raw_data.min_mz = std::numeric_limits<double>::infinity();
//     raw_data.max_mz = -std::numeric_limits<double>::infinity();
//     raw_data.min_rt = std::numeric_limits<double>::infinity();
//     raw_data.max_rt = -std::numeric_limits<double>::infinity();
//     raw_data.resolution_ms1 = resolution_ms1;
//     raw_data.resolution_msn = resolution_msn;
//     raw_data.reference_mz = reference_mz;
//     raw_data.fwhm_rt = 0;  // TODO(alex): Should this be passed as well?
//     raw_data.basicSpectra = {};
//     raw_data.retention_times = {};
//     // TODO(alex): Can we automatically detect the instrument type and set
//     // resolution from the header?

//     // Load the MzXML file
//     if (!mzfh.load(input_file.c_str())) {
//         std::cerr << "Error: Could not load file " << input_file << std::endl;
//         return std::nullopt;
//     }

//     // Measure the time to retrieve the SpectrumIndex
//     auto start_time = std::chrono::high_resolution_clock::now(); // Start timing

//     // Retrieve the SpectrumIndex
//     std::vector<cindex>* spectrumIndex = mzfh.getSpectrumIndex();

//     auto end_time = std::chrono::high_resolution_clock::now();   // End timing

//     // Calculate and print the elapsed time in milliseconds
//     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
//     std::cout << "Time taken to execute mzfh.getSpectrumIndex(): " << duration << " ms" << std::endl;

//     if (!spectrumIndex || spectrumIndex->empty()) {
//         std::cerr << "Warning: SpectrumIndex is empty or unavailable. Attempting to create it using scan range." << std::endl;

//         // Ensure lowScan and highScan are valid
//         int low_scan = mzfh.lowScan();
//         int high_scan = mzfh.highScan();

//         if (low_scan == -1 || high_scan == -1 || low_scan > high_scan) {
//             std::cerr << "Error: Invalid scan range. Unable to create SpectrumIndex." << std::endl;
//             return std::nullopt;
//         }

//         // Create SpectrumIndex based on the scan range
//         spectrumIndex = new std::vector<cindex>();
//         for (int scan_num = low_scan; scan_num <= high_scan; ++scan_num) {
//             cindex index;
//             index.scanNum = scan_num;
//             // Assuming an offset calculation can be derived; if not, set to 0 or handle appropriately
//             index.offset = 0; 
//             spectrumIndex->push_back(index);
//         }
//     }

//     // Iterate over the SpectrumIndex and read each spectrum

//     for (const auto& index : *spectrumIndex) {
//         if (mzfh.readSpectrumHeader(index.scanNum)) {
//             if (s.getMSLevel() == 1) { // Check the MS level from the header
//                 //ms1_indices.push_back(index.scanNum);
//                 mzfh.readSpectrum(index.scanNum);
//                 raw_data.basicSpectra.push_back(s);
                
//             }

// 		}
//     } // end reading spectrum loop
    
//     return raw_data;
// }

// Read an entire mzxml file into the RawData::RawData data structure filtering usin mstoolkit libraries1
// std::optional<RawData::RawMSDataS> XmlReader::read_msdatas(
//     std::string &input_file, double min_mz, double max_mz, double min_rt,
//     double max_rt, Instrument::Type instrument_type, double resolution_ms1,
//     double resolution_msn, double reference_mz, Polarity::Type polarity,
//     size_t ms_level) {

//     MSToolkit::Spectrum s;
//     // mz file handler
//     MSReader mzfh;
    
//     RawData::RawMSDataS raw_data = {};
//     raw_data.instrument_type = instrument_type;
//     raw_data.min_mz = std::numeric_limits<double>::infinity();
//     raw_data.max_mz = -std::numeric_limits<double>::infinity();
//     raw_data.min_rt = std::numeric_limits<double>::infinity();
//     raw_data.max_rt = -std::numeric_limits<double>::infinity();
//     raw_data.resolution_ms1 = resolution_ms1;
//     raw_data.resolution_msn = resolution_msn;
//     raw_data.reference_mz = reference_mz;
//     raw_data.fwhm_rt = 0;  // TODO(alex): Should this be passed as well?
//     raw_data.Spectra = {};
//     raw_data.retention_times = {};
//     // TODO(alex): Can we automatically detect the instrument type and set
//     // resolution from the header?

//     mzfh.addFilter(MS1);

//     // Iterate over the SpectrumIndex and read each spectrum
//     bool next = mzfh.readFile(input_file.c_str(),s);
//     while (next){
//         raw_data.Spectra.push_back(s);
//         next = mzfh.nextSpectrum(s);
//     }
  
//     return raw_data;
// }


// Read an entire mzxml file into the RawData::RawData data structure filtering using mstoolkit libraries
std::optional<RawData::RawData> XmlReader::read_msdata(
    std::string &input_file, double min_mz, double max_mz, double min_rt,
    double max_rt, Instrument::Type instrument_type, double resolution_ms1,
    double resolution_msn, double reference_mz, Polarity::Type polarity,
    size_t ms_level) {

    MSToolkit::Spectrum s;
    // mz file handler
    MSReader mzfh;
    
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

    if (ms_level < 1 || ms_level > 10) { 
        throw std::out_of_range("Invalid ms_level");
        return std::nullopt;
    }

    // Convert ms_level to the appropriate MSSpectrumType
    MSSpectrumType filter_type = static_cast<MSSpectrumType>(ms_level - 1);

    mzfh.addFilter(filter_type);

    // Iterate over the SpectrumIndex and read each spectrum
    bool next = mzfh.readFile(input_file.c_str(),s);

    if (!next) {
        std::cerr << "Error: Failed to read file: " << input_file << std::endl;
        // Optionally throw an exception or handle the error as needed
        throw std::runtime_error("Failed to read file: " + input_file);
    }
    
    while (next){
        auto scan = read_scan(s);
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
        next = mzfh.nextSpectrum(s);
    }
    
    return raw_data;
}

// Read an entire mzxml file into the RawData::RawData data structure filtering usin mstoolkit libraries1
// std::optional<RawData::RawData> XmlReader::read_mzxml2(
//     std::string &input_file, double min_mz, double max_mz, double min_rt,
//     double max_rt, Instrument::Type instrument_type, double resolution_ms1,
//     double resolution_msn, double reference_mz, Polarity::Type polarity,
//     size_t ms_level) {

//     BasicSpectrum s;
//     // mz file handler
//     MzParser mzfh(&s);
    
//     RawData::RawData raw_data = {};
//     raw_data.instrument_type = instrument_type;
//     raw_data.min_mz = std::numeric_limits<double>::infinity();
//     raw_data.max_mz = -std::numeric_limits<double>::infinity();
//     raw_data.min_rt = std::numeric_limits<double>::infinity();
//     raw_data.max_rt = -std::numeric_limits<double>::infinity();
//     raw_data.resolution_ms1 = resolution_ms1;
//     raw_data.resolution_msn = resolution_msn;
//     raw_data.reference_mz = reference_mz;
//     raw_data.fwhm_rt = 0;  // TODO(alex): Should this be passed as well?
//     raw_data.scans = {};
//     raw_data.retention_times = {};
//     // TODO(alex): Can we automatically detect the instrument type and set
//     // resolution from the header?

//     // Load the MzXML file
//     if (!mzfh.load(input_file.c_str())) {
//         std::cerr << "Error: Could not load file " << input_file << std::endl;
//         return std::nullopt;
//     }

//     // Retrieve the SpectrumIndex
//     std::vector<cindex>* spectrumIndex = mzfh.getSpectrumIndex();
//     if (!spectrumIndex || spectrumIndex->empty()) {
//         std::cerr << "Warning: SpectrumIndex is empty or unavailable. Attempting to create it using scan range." << std::endl;

//         // Ensure lowScan and highScan are valid
//         int low_scan = mzfh.lowScan();
//         int high_scan = mzfh.highScan();

//         if (low_scan == -1 || high_scan == -1 || low_scan > high_scan) {
//             std::cerr << "Error: Invalid scan range. Unable to create SpectrumIndex." << std::endl;
//             return std::nullopt;
//         }

//         // Create SpectrumIndex based on the scan range
//         spectrumIndex = new std::vector<cindex>();
//         for (int scan_num = low_scan; scan_num <= high_scan; ++scan_num) {
//             cindex index;
//             index.scanNum = scan_num;
//             // Assuming an offset calculation can be derived; if not, set to 0 or handle appropriately
//             index.offset = 0; 
//             spectrumIndex->push_back(index);
//         }
//     }

//     // Iterate over the SpectrumIndex and read each spectrum

//     for (const auto& index : *spectrumIndex) {
//         // Assuming cindex has meaningful fields (e.g., scan number, offset, etc.)
//         // std::cout << "Scan Number: " << index.scanNum << ", Offset: " << index.offset << std::endl;

// 		// if(index.scanNum<mzfh.lowScan() || index.scanNum>mzfh.highScan()) {
// 		// 	cout << "Bad number! BOOOOO!" << endl;
// 		// } else {
//     		if(!mzfh.readSpectrum(index.scanNum)) 
//                 cout << "Spectrum number not in file." << endl;
//     		else {
//                 auto scan = read_mzxml_scan(s);
//                 if (scan.num_points != 0) {
//                     raw_data.scans.push_back(scan);
//                     raw_data.retention_times.push_back(scan.retention_time);
//                     if (scan.retention_time < raw_data.min_rt)
//                         raw_data.min_rt = scan.retention_time;
//                     if (scan.retention_time > raw_data.max_rt) 
//                         raw_data.max_rt = scan.retention_time;
//                     if (scan.mz[0] < raw_data.min_mz)
//                         raw_data.min_mz = scan.mz[0];
//                     if (scan.mz[scan.mz.size() - 1] > raw_data.max_mz) 
//                         raw_data.max_mz = scan.mz[scan.mz.size() - 1];
//                 }
//             }    
// 		//}
//     } // end reading spectrum loop
    
//     return raw_data;
// }


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
        while (stream.good() && !stream.eof()) {
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
                while (stream.good() && !stream.eof()) {
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

                // Find whether or not the data is compressed.
                bool compressed = peak_attributes["compressionType"] == "zlib";

                // Extract the peaks from the data.
                auto data = XmlReader::read_data(stream);
                if (!data) {
                    return {};
                }

                // Decode base64-encoded string to raw data.
                std::vector<uint8_t> raw_data;
                Base64::decode_base64(data.value(), raw_data);

                if (compressed) {
                    // Calculate amount of bytes in decompressed data.
                    size_t decompressed_len = num_points * 2 * (precision / 8);
                    std::vector<uint8_t> decompressed_data;

                    // Decompress data.
                    int status = Compression::inflate(
                        raw_data, decompressed_data, decompressed_len);

                    // Check status after decompression.
                    if (status != Z_OK) {
                        return {};
                    }

                    raw_data = decompressed_data;
                }

                // Interpret raw data.
                double intensity_sum = 0;
                double max_intensity = 0;

                scan.mz.resize(num_points);
                scan.intensity.resize(num_points);
                size_t scan_size = 0;

                size_t offset = 0;
                for (size_t i = 0; i < num_points; ++i) {
                    double mz = 0;
                    double intensity = 0;
                    // Interpret mz and intensity values.
                    if (precision == 32) {
                        mz = Base64::interpret_float(raw_data, offset,
                                                     little_endian);
                        offset += 4;
                        intensity = Base64::interpret_float(raw_data, offset,
                                                            little_endian);
                        offset += 4;
                    } else if (precision == 64) {
                        mz = Base64::interpret_double(raw_data, offset,
                                                      little_endian);
                        offset += 8;
                        intensity = Base64::interpret_double(raw_data, offset,
                                                             little_endian);
                        offset += 8;
                    }

                    // We don't need to extract the peaks when we are not
                    // inside the mz bounds or contain no value.
                    if (mz < min_mz || mz > max_mz || intensity == 0) {
                        continue;
                    }

                    if (intensity > max_intensity) {
                        max_intensity = intensity;
                    }
                    intensity_sum += intensity;

                    scan.mz[scan_size] = mz;
                    scan.intensity[scan_size] = intensity;
                    ++scan_size;
                }
                // Resize to number of elements that are actually included.
                scan.mz.resize(scan_size);
                scan.intensity.resize(scan_size);

                // Shrink capacity
                scan.mz.shrink_to_fit();
                scan.intensity.shrink_to_fit();

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
        while (stream.good() && next_tag) {
            if (next_tag.value().name == "scan" && next_tag.value().closed) {
                break;
            }
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
    while (stream.good() && !stream.eof()) {
        auto tag = XmlReader::read_tag(stream);
        if (!tag) {
            continue;
        }
        auto atts = tag.value().attributes;

        if (tag.value().name == "msRun" && tag.value().closed) {
            break;
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

std::optional<RawData::RawData> XmlReader::_read_mzml(
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
    while (stream.good() && !stream.eof()) {
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
            while (stream.good() && !stream.eof()) {
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
                    auto precursor_attributes = tag.value().attributes;
                    std::string spectrumRef =
                        precursor_attributes["spectrumRef"];
                    // Find scan number.
                    size_t scanIdx = spectrumRef.find("scan=") + 5;
                    scan.precursor_information.scan_number =
                        std::stoull(spectrumRef.substr(scanIdx));

                    scan.precursor_information.charge = 0;
                    scan.precursor_information.mz = 0.0;
                    scan.precursor_information.window_wideness = 0.0;
                    scan.precursor_information.intensity = 0.0;
                    scan.precursor_information.activation_method =
                        ActivationMethod::UNKNOWN;
                    while (stream.good() && !stream.eof()) {
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
                    // Uncompressed: false, Zlib compression: true.
                    bool compressed = false;
                    // mz: 0, intensity: 1
                    int type = -1;
                    std::optional<std::string> data;
                    while (stream.good() && !stream.eof()) {
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
                            // Compression.
                            if (accession == "MS:1000574") {
                                compressed = true;
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
                        // decode data.
                        std::vector<uint8_t> binary_data;
                        Base64::decode_base64(data.value(), binary_data);
                        if (compressed) {
                            std::vector<uint8_t> decompressed_data;

                            // Decompress data, set decompressed length to 0
                            // (unknown).
                            int status = Compression::inflate(
                                binary_data, decompressed_data, 0);

                            // Check status after decompression.
                            if (status != Z_OK) {
                                return raw_data;
                            }

                            binary_data = decompressed_data;
                        }

                        size_t num_points =
                            binary_data.size() / (precision / 8);

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

                        // offset for interpretation.
                        size_t offset = 0;
                        for (size_t i = 0; i < num_points; ++i) {
                            double value = 0;
                            if (precision == 32) {
                                value = Base64::interpret_float(binary_data,
                                                                offset, true);
                                offset += 4;
                            } else if (precision == 64) {
                                value = Base64::interpret_double(binary_data,
                                                                 offset, true);
                                offset += 8;
                            }
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
    auto tag = std::optional<Tag>(Tag{});
    if (!stream.good() || stream.eof()) {
        return std::nullopt;
    }

    // Store the tag contents in a buffer for further processing.
    std::string buffer;

    while (stream.good() && !stream.eof() && stream.get() != '<')
        ;
    std::getline(stream, buffer, '>');

    if (buffer.empty()) {
        return std::nullopt;
    }

    // Check if this is a closing tag for a previous one.
    if (buffer[0] == '/') {
        tag->closed = true;
        // Read tag name
        size_t end_name = buffer.size();
        for (size_t i = 1; i < buffer.size(); ++i) {
            if (std::isspace(buffer[i])) {
                end_name = i;
                break;
            }
        }
        tag->name = buffer.substr(1, end_name);
        return tag;
    }

    // Read tag name
    size_t start_ptr = buffer.size();
    size_t end_ptr = buffer.size();
    for (size_t i = 1; i < buffer.size(); ++i) {
        if (std::isspace(buffer[i])) {
            start_ptr = i;
            break;
        }
    }
    tag->name = buffer.substr(0, start_ptr);

    // Check if this is a self-closing tag.
    if (buffer[buffer.size() - 1] == '/') {
        end_ptr = buffer.size() - 1;
        tag->closed = true;
    }

    // Read attributes
    for (size_t i = start_ptr; i < end_ptr; ++i) {
        // Find attribute name.
        if (std::isspace(buffer[i])) {
            continue;
        }
        size_t ptr = i;
        while (buffer[ptr] != '=') {
            // Malformed xml.
            if (ptr == end_ptr) {
                break;
            }
            ptr++;
        }
        if (ptr == end_ptr) {
            break;
        }
        std::string att_name = buffer.substr(i, ptr - i);

        // Find attribute value.
        i = ptr;
        ptr--;
        bool beg = false;
        bool end = false;
        bool escaped = false;
        while (!(beg && end)) {
            // Malformed xml.
            if (ptr == end_ptr) {
                return std::nullopt;
            }
            if (escaped) {
                escaped = false;
                ptr++;
                continue;
            }

            char c = buffer[ptr];
            if (c == '\\') {
                escaped = true;
            }
            if (c == '"') {
                if (!beg) {
                    beg = true;
                    ptr++;
                    continue;
                }
                if (!end) {
                    end = true;
                    break;
                }
            }
            ptr++;
        }
        std::string att_value = buffer.substr(i + 2, ptr - i - 2);
        i = ptr;

        tag->attributes[att_name] = att_value;
    }

    return tag;
}

IdentData::IdentData XmlReader::read_mzidentml(std::istream &stream,
                                               bool ignore_decoy,
                                               bool require_threshold,
                                               bool max_rank_only,
                                               double min_mz, double max_mz,
                                               double min_rt, double max_rt) {
    IdentData::IdentData ident_data = {};

    // Find the DBSequences, Peptides and PeptideEvidence in the
    // SequenceCollection tag.
    while (stream.good() && !stream.eof()) {
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
                while (stream.good() && !stream.eof()) {
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
            while (stream.good() && !stream.eof()) {
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
                    while (stream.good() && !stream.eof()) {
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
    while (stream.good() && !stream.eof()) {
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
        while (stream.good() && !stream.eof()) {
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

                spectrum_matches.push_back(spectrum_match);
            }
        }

        if (max_rank_only) {
            if (spectrum_matches.empty()) {
                continue;
            }
            IdentData::SpectrumMatch selected_spectrum = {};
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
            if (selected_spectrum.experimental_mz >= min_mz &&
                selected_spectrum.experimental_mz <= max_mz &&
                selected_spectrum.retention_time >= min_rt &&
                selected_spectrum.retention_time <= max_rt) {
                ident_data.spectrum_matches.push_back(selected_spectrum);
            }
        } else {
            // Update retention time on the provisional spectrum_matches list
            // and push each element to the list of PSM.
            for (auto &spectrum_match : spectrum_matches) {
                spectrum_match.retention_time = retention_time;
                if (spectrum_match.experimental_mz >= min_mz &&
                    spectrum_match.experimental_mz <= max_mz &&
                    spectrum_match.retention_time >= min_rt &&
                    spectrum_match.retention_time <= max_rt) {
                    ident_data.spectrum_matches.push_back(spectrum_match);
                }
            }
        }
    }
    return ident_data;
}
