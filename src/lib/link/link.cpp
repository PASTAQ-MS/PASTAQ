#include <algorithm>

#include "link/link.hpp"
#include "utils/search.hpp"

std::vector<Link::LinkedMsms> Link::link_peaks(
        const std::vector<Centroid::Peak> &peaks,
        const RawData::RawData &raw_data, double n_sig_mz, double n_sig_rt) {
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
        double min_mz = event_mz - n_sig_mz * theoretical_sigma_mz;
        double max_mz = event_mz + n_sig_mz * theoretical_sigma_mz;
        size_t min_j = Search::lower_bound(index_mzs, min_mz);
        double min_distance = std::numeric_limits<double>::infinity();
        size_t peak_id = 0;
        for (size_t j = min_j; j < index_ids.size(); ++j) {
            size_t i = index_ids[j];
            if (peaks[i].fitted_mz > max_mz) {
                break;
            }
            double a = (event_mz - peaks[i].fitted_mz) / peaks[i].fitted_sigma_mz;
            double b = (event_rt - peaks[i].fitted_rt) / peaks[i].fitted_sigma_rt;
            double distance = std::sqrt(a * a + b * b);
            if (distance < min_distance) {
                min_distance = distance;
                peak_id = i;
            }
        }

        // Check if linked event is within n sigma of the minimum distance
        // peak.
        double roi_min_mz =
            peaks[peak_id].fitted_mz - n_sig_mz * peaks[peak_id].fitted_sigma_mz;
        double roi_max_mz =
            peaks[peak_id].fitted_mz + n_sig_mz * peaks[peak_id].fitted_sigma_mz;
        double roi_min_rt =
            peaks[peak_id].fitted_rt - n_sig_rt * peaks[peak_id].fitted_sigma_rt;
        double roi_max_rt =
            peaks[peak_id].fitted_rt + n_sig_rt * peaks[peak_id].fitted_sigma_rt;
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
        const IdentData::IdentData &ident_data, 
        const RawData::RawData &raw_data, double n_sig_mz, double n_sig_rt) {
    // Index the msms list by m/z.
    struct MsmsIndex {
        uint64_t id;
        double mz;
    };
    auto indices = std::vector<MsmsIndex>(raw_data.scans.size());
    for (size_t i = 0; i < raw_data.scans.size(); ++i) {
        const auto &scan = raw_data.scans[i];
        if (scan.ms_level != 2) {
            continue;
        }
        indices[i] = {i, scan.precursor_information.mz};
    }

    // Sort mz index by mz.
    std::sort(indices.begin(), indices.end(),
              [](const MsmsIndex &a, const MsmsIndex &b) -> bool {
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

    // Perform linkage of ms/ms events to closest PSM.
    std::vector<Link::LinkedMsms> link_table;
    for (size_t k = 0; k < ident_data.spectrum_matches.size(); ++k) {
        const auto &psm = ident_data.spectrum_matches[k];

        size_t psm_index = k;
        double psm_mz = psm.experimental_mz;
        double psm_rt = psm.retention_time;
        double theoretical_sigma_mz = RawData::fwhm_to_sigma(
            RawData::theoretical_fwhm(raw_data, psm_mz));
        double theoretical_sigma_rt = RawData::fwhm_to_sigma(raw_data.fwhm_rt);

        // Find min_mz and loop until we reach the max_mz.
        double min_mz = psm_mz - n_sig_mz * theoretical_sigma_mz;
        double max_mz = psm_mz + n_sig_mz * theoretical_sigma_mz;
        size_t min_j = Search::lower_bound(index_mzs, min_mz);
        double min_distance = std::numeric_limits<double>::infinity();
        size_t selected_i = 0;
        for (size_t j = min_j; j < index_ids.size(); ++j) {
            size_t i = index_ids[j];
            const auto &scan = raw_data.scans[i];
            double scan_mz = scan.precursor_information.mz;
            double scan_rt= scan.retention_time;

            if (scan.precursor_information.mz > max_mz) {
                break;
            }
            double a = (psm_mz - scan_mz) / scan_mz;
            double b = (psm_rt - scan_rt) / scan_rt;
            double distance = std::sqrt(a * a + b * b);
            if (distance < min_distance) {
                min_distance = distance;
                selected_i = i;
            }
        }

        if (min_distance == std::numeric_limits<double>::infinity()) {
            continue;
        }

        const auto &scan = raw_data.scans[selected_i];
        double scan_mz = scan.precursor_information.mz;
        double scan_rt= scan.retention_time;

        // Check if linked event is within 3 sigma of the minimum distance scan.
        double roi_min_mz = scan_mz - n_sig_mz * theoretical_sigma_mz;
        double roi_max_mz = scan_mz + n_sig_mz * theoretical_sigma_mz;
        double roi_min_rt = scan_rt - n_sig_rt * theoretical_sigma_rt;
        double roi_max_rt = scan_rt + n_sig_rt * theoretical_sigma_rt;
        if (psm_mz < roi_min_mz || psm_mz > roi_max_mz ||
            psm_rt < roi_min_rt || psm_rt > roi_max_rt) {
            continue;
        }

        link_table.push_back({psm_index, scan.scan_number, selected_i, min_distance});
    }

    // Sort link_table by msms_id.
    std::sort(
        link_table.begin(), link_table.end(),
        [](const Link::LinkedMsms &a, const Link::LinkedMsms &b) -> bool {
            return (a.msms_id < b.msms_id);
        });
    return link_table;
}

// Experimental PSM to peak linkage based on theoretical m/z instead of
// experimental m/z or MSMS id (scan index).
std::vector<Link::LinkedPsm> Link::link_psm(
        const IdentData::IdentData &ident_data,
        const std::vector<Centroid::Peak> &peaks,
        const RawData::RawData &raw_data, double n_sig_mz, double n_sig_rt) {
    // Index the peak list by m/z.
    struct PeakIndex {
        uint64_t id;
        double mz;
    };
    auto indices = std::vector<PeakIndex>(peaks.size());
    for (size_t i = 0; i < peaks.size(); ++i) {
        indices[i] = {i, peaks[i].fitted_mz};
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
    std::vector<Link::LinkedPsm> link_table;
    for (size_t k = 0; k < ident_data.spectrum_matches.size(); ++k) {
        const auto &psm = ident_data.spectrum_matches[k];

        size_t psm_index = k;
        double psm_mz = psm.theoretical_mz;
        double psm_rt = psm.retention_time;
        // FIXME: This is ugly, we don't need rawdata anywere else in this
        // function. We should decouple theoretical mz calculation from the raw
        // data...
        double theoretical_sigma_mz = RawData::fwhm_to_sigma(
            RawData::theoretical_fwhm(raw_data, psm_mz));

        // Find min_mz and loop until we reach the max_mz.
        double min_mz = psm_mz - n_sig_mz * theoretical_sigma_mz;
        double max_mz = psm_mz + n_sig_mz * theoretical_sigma_mz;
        size_t min_j = Search::lower_bound(index_mzs, min_mz);
        double min_distance = std::numeric_limits<double>::infinity();
        size_t selected_i = 0;
        for (size_t j = min_j; j < index_ids.size(); ++j) {
            size_t i = index_ids[j];
            if (peaks[i].fitted_mz > max_mz) {
                break;
            }
            double a = (psm_mz - peaks[i].fitted_mz) / peaks[i].fitted_sigma_mz;
            double b = (psm_rt - peaks[i].fitted_rt) / peaks[i].fitted_sigma_rt;
            double distance = std::sqrt(a * a + b * b);
            if (distance < min_distance) {
                min_distance = distance;
                selected_i = i;
            }
        }
        const auto &peak = peaks[selected_i];

        // Check if linked event is within 3 sigma of the minimum distance
        // peak.
        double roi_min_mz = peak.fitted_mz - n_sig_mz * peak.fitted_sigma_mz;
        double roi_max_mz = peak.fitted_mz + n_sig_mz * peak.fitted_sigma_mz;
        double roi_min_rt = peak.fitted_rt - n_sig_rt * peak.fitted_sigma_rt;
        double roi_max_rt = peak.fitted_rt + n_sig_rt * peak.fitted_sigma_rt;
        if (psm_mz < roi_min_mz || psm_mz > roi_max_mz ||
            psm_rt < roi_min_rt || psm_rt > roi_max_rt) {
            continue;
        }

        if (min_distance != std::numeric_limits<double>::infinity()) {
            link_table.push_back({peak.id, psm_index, min_distance});
        }
    }

    // Sort link_table by peak_id.
    std::sort(
        link_table.begin(), link_table.end(),
        [](const Link::LinkedPsm &a, const Link::LinkedPsm &b) -> bool {
            return (a.peak_id < b.peak_id) ||
                   ((a.peak_id == b.peak_id) && (a.distance < b.distance));
        });
    return link_table;
}
