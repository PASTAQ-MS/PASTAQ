#include <algorithm>

#include "link/link.hpp"
#include "utils/search.hpp"

std::vector<Link::LinkedMsms> Link::link_peaks(
    const std::vector<Centroid::Peak> &peaks,
    const RawData::RawData &raw_data) {
    // Index the peak list by m/z.
    struct PeakIndex {
        uint64_t id;
        double mz;
    };
    auto indices = std::vector<PeakIndex>(peaks.size());
    for (size_t i = 0; i < peaks.size(); ++i) {
        indices[i] = {peaks[i].id, peaks[i].fitted_mz};
    }

    // Sort mz index by mz.
    std::sort(indices.begin(), indices.end(),
              [](const PeakIndex &a, const PeakIndex &b) -> bool {
                  return a.mz < b.mz;
              });

    // Flatten index into separate arrays.
    auto index_ids = std::vector<size_t>(indices.size());
    auto index_mzs = std::vector<double>(indices.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        index_ids[i] = indices[i].id;
        index_mzs[i] = indices[i].mz;
    }
    indices.clear();

    // Perform linkage of ms/ms events to closest peak.
    std::vector<Link::LinkedMsms> link_table;
    for (size_t k = 0; k < raw_data.scans.size(); ++k) {
        const auto &scan = raw_data.scans[k];
        if (scan.ms_level != 2) {
            continue;
        }
        size_t event_id = scan.scan_number;
        double event_mz = scan.precursor_information.mz;
        double event_rt = scan.retention_time;
        double theoretical_sigma_mz = RawData::fwhm_to_sigma(
            RawData::theoretical_fwhm(raw_data, event_mz));

        // Find min_mz and loop until we reach the max_mz.
        double min_mz = event_mz - 3 * theoretical_sigma_mz;
        double max_mz = event_mz + 3 * theoretical_sigma_mz;
        size_t min_j = Search::lower_bound(index_mzs, min_mz);
        double min_distance = std::numeric_limits<double>::infinity();
        size_t peak_id = 0;
        for (size_t j = min_j; j < index_ids.size(); ++j) {
            size_t i = index_ids[j];
            if (peaks[i].fitted_mz > max_mz) {
                break;
            }
            double a = event_mz - peaks[i].fitted_mz;
            double b = event_rt - peaks[i].fitted_rt;
            double distance = std::sqrt(a * a + b * b);
            if (distance < min_distance) {
                min_distance = distance;
                peak_id = i;
            }
        }

        // Check if linked event is within 10 sigma of the minimum distance
        // peak.
        double roi_min_mz =
            peaks[peak_id].fitted_mz - 10 * peaks[peak_id].fitted_sigma_mz;
        double roi_max_mz =
            peaks[peak_id].fitted_mz + 10 * peaks[peak_id].fitted_sigma_mz;
        double roi_min_rt =
            peaks[peak_id].fitted_rt - 10 * peaks[peak_id].fitted_sigma_rt;
        double roi_max_rt =
            peaks[peak_id].fitted_rt + 10 * peaks[peak_id].fitted_sigma_rt;
        if (event_mz < roi_min_mz || event_mz > roi_max_mz ||
            event_rt < roi_min_rt || event_rt > roi_max_rt) {
            continue;
        }

        if (min_distance != std::numeric_limits<double>::infinity()) {
            link_table.push_back({peak_id, event_id, k, min_distance});
        }
    }

    // Sort link_table by entity_id.
    std::sort(
        link_table.begin(), link_table.end(),
        [](const Link::LinkedMsms &a, const Link::LinkedMsms &b) -> bool {
            return (a.entity_id < b.entity_id) ||
                   ((a.entity_id == b.entity_id) && (a.distance < b.distance));
        });
    return link_table;
}

std::vector<Link::LinkedMsms> Link::link_idents(
    const IdentData::IdentData &ident_data, const RawData::RawData &raw_data) {
    // Index the spectrumid list by m/z.
    struct SpectrumIdIndex {
        uint64_t id;
        double mz;
    };
    auto &spectrum_ids = ident_data.spectrum_ids;
    if (spectrum_ids.empty()) {
        return {};
    }
    auto indices = std::vector<SpectrumIdIndex>(spectrum_ids.size());
    for (size_t i = 0; i < spectrum_ids.size(); ++i) {
        indices[i] = {i, spectrum_ids[i].experimental_mz};
    }

    // Sort mz index by mz.
    std::sort(indices.begin(), indices.end(),
              [](const SpectrumIdIndex &a, const SpectrumIdIndex &b) -> bool {
                  return a.mz < b.mz;
              });

    // Flatten index into separate arrays.
    auto index_ids = std::vector<size_t>(indices.size());
    auto index_mzs = std::vector<double>(indices.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        index_ids[i] = indices[i].id;
        index_mzs[i] = indices[i].mz;
    }
    indices.clear();

    // Perform linkage of ms/ms events to closest spectrum_id.
    std::vector<Link::LinkedMsms> link_table;
    for (size_t k = 0; k < raw_data.scans.size(); ++k) {
        const auto &scan = raw_data.scans[k];
        if (scan.ms_level != 2) {
            continue;
        }
        size_t event_id = scan.scan_number;
        double event_mz = scan.precursor_information.mz;
        double event_rt = scan.retention_time;
        double theoretical_sigma_mz = RawData::fwhm_to_sigma(
            RawData::theoretical_fwhm(raw_data, event_mz));

        // Find min_mz and loop until we reach the max_mz.
        double min_mz = event_mz - 1 * theoretical_sigma_mz;
        double max_mz = event_mz + 1 * theoretical_sigma_mz;
        size_t min_j = Search::lower_bound(index_mzs, min_mz);
        double min_distance = std::numeric_limits<double>::infinity();
        size_t spectrum_id_id = 0;
        for (size_t j = min_j; j < index_ids.size(); ++j) {
            size_t i = index_ids[j];
            if (spectrum_ids[i].experimental_mz > max_mz) {
                break;
            }
            double a = event_mz - spectrum_ids[i].experimental_mz;
            double b = event_rt - spectrum_ids[i].retention_time;
            double distance = std::sqrt(a * a + b * b);
            if (distance < min_distance) {
                min_distance = distance;
                spectrum_id_id = i;
            }
        }

        if (min_distance != std::numeric_limits<double>::infinity()) {
            link_table.push_back({spectrum_id_id, event_id, k, min_distance});
        }
    }

    // Sort link_table by entity_id.
    std::sort(
        link_table.begin(), link_table.end(),
        [](const Link::LinkedMsms &a, const Link::LinkedMsms &b) -> bool {
            return (a.entity_id < b.entity_id) ||
                   ((a.entity_id == b.entity_id) && (a.distance < b.distance));
        });
    return link_table;
}
