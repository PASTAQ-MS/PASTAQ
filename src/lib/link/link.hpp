#ifndef LINK_LINKMSMS_HPP
#define LINK_LINKMSMS_HPP

#include <iostream>
#include <vector>

#include "centroid/centroid.hpp"
#include "raw_data/raw_data.hpp"

namespace Link {
// This table can store linkage for both peak-msms and spectrum_id-msms. In
// both cases, entity_id will store the index for the peak or spectrum_id lists
// respectively.
// TODO(alex): It is convenient for being able to return the same object, but
// since they have slightly different meanings it could be a bit confusing.
struct LinkedMsms {
    uint64_t entity_id;
    uint64_t msms_id;
    uint64_t scan_index;
    double distance;
};
// TODO(alex): This needs more documentation

// NOTE: This algorithm relies on the peak vector to be sorted by id/height.
std::vector<LinkedMsms> link_peaks(const std::vector<Centroid::Peak> &peaks,
        const RawData::RawData &raw_data, double n_sig_mz, double n_sig_rt);
std::vector<LinkedMsms> link_idents(const IdentData::IdentData &ident_data,
        const RawData::RawData &raw_data, double n_sig_mz, double n_sig_rt);

struct LinkedPsm {
    uint64_t peak_id;
    uint64_t psm_index;
    double distance;
};
std::vector<LinkedPsm> link_psm(const IdentData::IdentData &ident_data,
        const std::vector<Centroid::Peak> &peaks,
        const RawData::RawData &raw_data, double n_sig_mz, double n_sig_rt);

}  // namespace Link

#endif /* LINK_LINKMSMS_HPP */
