#include "MIDAs/MIDAs.h"

#include "feature_detection/feature_detection.hpp"

std::tuple<std::vector<double>, std::vector<double>>
FeatureDetection::theoretical_isotopes_peptide(std::string sequence,
                                               int8_t charge_state,
                                               double min_perc) {
    auto midas = MIDAs(charge_state = charge_state);
    // NOTE: C00: Unmodified cysteine.
    //       C31: Carboxymethylation or Iodoacetic acid.
    //       C32: Carbamidomethylation or Iodoacetamide.
    //       C33: Pyridylethylation.
    midas.Initialize_Elemental_Composition(sequence, "C00", "H", "OH", 1);
    auto isotopes = midas.Coarse_Grained_Isotopic_Distribution();
    auto full_mzs = std::vector<double>(isotopes.size());
    auto probs = std::vector<double>(isotopes.size());
    double max_prob = 0;
    for (size_t i = 0; i < isotopes.size(); ++i) {
        full_mzs[i] = isotopes[i].mw / charge_state;
        probs[i] = isotopes[i].prob;
        if (probs[i] > max_prob) {
            max_prob = probs[i];
        }
    }
    // Normalize the probabilities to 0-1.
    std::vector<double> mzs;
    std::vector<double> perc;
    for (size_t i = 0; i < isotopes.size(); ++i) {
        probs[i] = probs[i] / max_prob;
        if (probs[i] > min_perc) {
            mzs.push_back(full_mzs[i]);
            perc.push_back(probs[i]);
        }
    }
    return {mzs, perc};
}

std::optional<FeatureDetection::Feature> FeatureDetection::build_feature(
    const std::vector<bool> &peaks_in_use,
    const std::vector<Centroid::Peak> &peaks,
    const std::vector<Search::KeySort<double>> &peaks_rt_key,
    const std::vector<double> &mzs, const std::vector<double> &perc,
    double tolerance_rt, double retention_time, double discrepancy_threshold) {
    // Basic sanitation.
    if (mzs.size() != perc.size() || mzs.size() == 0) {
        return std::nullopt;
    }

    // Find the peaks in range for matching.
    double min_rt = retention_time - tolerance_rt;
    double max_rt = retention_time + tolerance_rt;
    size_t min_j = Search::lower_bound(peaks_rt_key, min_rt);
    size_t max_j = peaks_rt_key.size();
    std::vector<Centroid::Peak> peaks_in_range;
    for (size_t j = min_j; j < max_j; ++j) {
        if (peaks_rt_key[j].sorting_key > max_rt) {
            break;
        }
        auto &peak = peaks[peaks_rt_key[j].index];
        if (peaks_in_use[peak.id]) {
            continue;
        }
        if ((peak.local_max_mz + peak.raw_roi_sigma_mz * 2) < mzs[0] ||
            (peak.local_max_mz - peak.raw_roi_sigma_mz * 2) >
                mzs[mzs.size() - 1]) {
            continue;
        }
        peaks_in_range.push_back(peak);
    }
    if (peaks_in_range.empty()) {
        return std::nullopt;
    }
    // Sort the peaks in range by mz for a faster search.
    std::sort(peaks_in_range.begin(), peaks_in_range.end(),
                     [](const Centroid::Peak &p1, const Centroid::Peak &p2) {
                         return (p1.local_max_mz < p2.local_max_mz);
                     });

    // Find the reference node and the list of candidates for each node.
    size_t reference_node_index = 0;
    bool found_reference_node = false;
    std::vector<std::vector<const Centroid::Peak *>> candidate_list;
    for (size_t k = 0; k < mzs.size(); ++k) {
        if (perc[k] == 1) {
            reference_node_index = k;
            found_reference_node = true;
        }
        double theoretical_mz = mzs[k];
        std::vector<const Centroid::Peak *> candidates_node;
        for (const auto &peak : peaks_in_range) {
            if ((peak.local_max_mz + peak.raw_roi_sigma_mz) > theoretical_mz &&
                (peak.local_max_mz - peak.raw_roi_sigma_mz) < theoretical_mz) {
                candidates_node.push_back(&peak);
            }
        }
        candidate_list.push_back(candidates_node);
    }
    if (!found_reference_node || candidate_list.empty()) {
        return std::nullopt;
    }

    // In case more than one peak is linked to the reference mz isotope,
    // the sequence with the less matching error should be selected. In
    // order to do so, the heights for each candidate must be normalized
    // by the reference isotope height.
    std::vector<const Centroid::Peak *> selected_candidates;
    std::vector<double> selected_candidates_norm_height;
    double min_total_distance = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < candidate_list[reference_node_index].size(); ++i) {
        const auto ref_candidate = candidate_list[reference_node_index][i];
        std::vector<const Centroid::Peak *> selected_candidates_for_reference;
        std::vector<double> selected_candidates_for_reference_norm_height;
        double total_distance = 0.0;
        // We need to do a forwards and backwards pass in order for
        // the algorithm to work, since we are be stopping the pattern matching
        // once we find a peak where the minimum height difference is greater
        // than the given discrepancy_threshold.
        //     reference_node_index->0           // Backwards
        //     reference_node_index->mzs.size()  // Forward
        // The peaks are reversed after the backwards pass, so that we keep the
        // proper ascending mz ordering.
        auto select_candidates = [&](size_t k) -> bool {
            if (k == reference_node_index) {
                selected_candidates_for_reference.push_back(ref_candidate);
                selected_candidates_for_reference_norm_height.push_back(1.0);
                return true;
            }

            // Find the best matching candidate for the selected reference. We
            // define here the best matching candidate as the candidate peak
            // with the minimum euclidean distance from the expected
            // theoretical peak.
            double theoretical_mz = mzs[k];
            double theoretical_percentage = perc[k];
            double min_distance = std::numeric_limits<double>::infinity();
            double selected_normalized_height = 0.0;
            const Centroid::Peak *selected_candidate;
            double selected_normalized_height_diff = 0.0;
            for (size_t j = 0; j < candidate_list[k].size(); ++j) {
                const auto candidate = candidate_list[k][j];
                double normalized_height = candidate->local_max_height /
                                           ref_candidate->local_max_height;
                double a = candidate->local_max_mz - theoretical_mz;
                double b = candidate->local_max_rt - retention_time;
                double c = normalized_height - theoretical_percentage;
                double distance = std::sqrt(a * a + b * b + c * c);
                if (distance < min_distance) {
                    min_distance = distance;
                    selected_normalized_height = normalized_height;
                    selected_candidate = candidate;
                    selected_normalized_height_diff = std::abs(c);
                }
            }

            if (selected_normalized_height_diff > discrepancy_threshold ||
                selected_normalized_height == 0.0) {
                return false;
            }
            selected_candidates_for_reference.push_back(selected_candidate);
            selected_candidates_for_reference_norm_height.push_back(
                selected_normalized_height);
            total_distance += min_distance;
            return true;
        };
        // Backwards from reference_node_index.
        for (size_t k = reference_node_index; k > 0; --k) {
            int success = select_candidates(k);
            if (success) {
                continue;
            }
            if (!success) {
                break;
            }
        }
        // Candidate at exactly 0.
        select_candidates(0);
        // Reverse list of candidates.
        std::reverse(selected_candidates_for_reference.begin(),
                     selected_candidates_for_reference.end());
        std::reverse(selected_candidates_for_reference_norm_height.begin(),
                     selected_candidates_for_reference_norm_height.end());
        // Forward from reference_node_index.
        for (size_t k = reference_node_index + 1; k < mzs.size(); ++k) {
            int success = select_candidates(k);
            if (success) {
                continue;
            }
            if (!success) {
                break;
            }
        }
        if (total_distance < min_total_distance) {
            selected_candidates = selected_candidates_for_reference;
            selected_candidates_norm_height =
                selected_candidates_for_reference_norm_height;
            min_total_distance = total_distance;
        }
    }
    if (selected_candidates.empty()) {
        return std::nullopt;
    }

    // Build the actual feature data.
    Feature feature = {};
    feature.rt = retention_time;
    // FIXME: Currently assuming that the minimum detected isotope is the
    // monoisotopic peak, but THIS MIGHT NOT BE THE CASE. For simplicity and to
    // keep the flow going I'll leave this for now, but must go back and FIX
    // it.
    feature.monoisotopic_mz = selected_candidates[0]->local_max_mz;
    feature.monoisotopic_height = selected_candidates[0]->local_max_height;
    // Find the weighted average mz and the average height of the selected
    // isotopes.
    feature.average_mz = 0.0;
    feature.total_height = 0.0;
    for (size_t i = 0; i < selected_candidates.size(); ++i) {
        auto candidate = selected_candidates[i];
        feature.total_height += candidate->local_max_height;
        feature.average_mz +=
            candidate->local_max_height * candidate->local_max_mz;
        feature.peak_ids.push_back(candidate->id);
    }
    if (feature.total_height == 0) {
        return std::nullopt;
    }
    feature.average_mz /= feature.total_height;
    feature.average_mz = feature.average_mz;
    return feature;
}

std::vector<FeatureDetection::Feature> FeatureDetection::feature_detection(
    const std::vector<Centroid::Peak> &peaks,
    const RawData::RawData &raw_data_ms2,
    const IdentData::IdentData &ident_data,
    const std::vector<Link::LinkedMsms> &link_table_msms,
    const std::vector<Link::LinkedMsms> &link_table_idents,
    double discrepancy_threshold) {
    std::vector<Feature> features;

    // The proposed algorithm goes as follows:
    // [X] 0.- Copy and sort the necessary vectors for the use of binary search
    //         (link_table_idents is sorted by msms, as that is the key being
    //         used for searching). The peaks array and link_table_msms array
    //         should have already been sorted by peak_id, which is what we
    //         want, as we are going to start matching peaks from highest
    //         intensity to lowest.
    // [X] 1.- For each linked peak on the link_table_msms, find it's associated
    //         entry on the link_table_idents.
    // [X] 2.- If the entry is found, use it to generate a theoretical isotopic
    //         distribution, otherwise, averagine will be generated.
    // [ ] 3.- We try to find the proposed peaks from the theoretical
    //         distribution on the peaks array (Maximum likelihood).
    // [ ] 4.- The peaks are marked as non available for future use. This means
    //         that this is a greedy algorithm.

    // Copy and sort key vectors.
    auto idents_msms_key =
        std::vector<Search::KeySort<uint64_t>>(link_table_idents.size());
    for (size_t i = 0; i < link_table_idents.size(); ++i) {
        idents_msms_key[i] = {i, link_table_idents[i].msms_id};
    }
    {
        auto sorting_key_func = [](const Search::KeySort<uint64_t> &p1,
                                   const Search::KeySort<uint64_t> &p2) {
            return (p1.sorting_key < p2.sorting_key);
        };
        std::sort(idents_msms_key.begin(), idents_msms_key.end(),
                         sorting_key_func);
    }
    auto peaks_rt_key = std::vector<Search::KeySort<double>>(peaks.size());
    for (size_t i = 0; i < peaks.size(); ++i) {
        peaks_rt_key[i] = {i, peaks[i].local_max_rt};
    }
    {
        auto sorting_key_func = [](const Search::KeySort<double> &p1,
                                   const Search::KeySort<double> &p2) {
            return (p1.sorting_key < p2.sorting_key);
        };
        std::sort(peaks_rt_key.begin(), peaks_rt_key.end(),
                         sorting_key_func);
    }
    // We use this variable to keep track of the peaks we have already linked.
    auto peaks_in_use = std::vector<bool>(peaks.size());

    // TODO: We should probably prioritize the MSMS events that HAVE an
    // identification, instead of just being intensity based only.
    // TODO: We are linking msms events independently, but we know that
    // multiple msms events can be linked to any given peak. If the selected
    // identification is the same there is no problem, and if there is
    // a majority of identifications they should result in a consensus.
    // However, if the identification is ambiguous, we can choose to use
    // averagine instead for that peak. To avoid duplication of effor, we
    // should resolve conflicts and perform feature matching based only on the
    // list of consensus candidates.
    for (const auto &linked_msms : link_table_msms) {
        // Find lower bound on the idents_msms_key.
        // auto i = lower_bound(idents_msms_key, 99999);
        auto i = lower_bound(idents_msms_key, linked_msms.msms_id);
        auto charge_state = raw_data_ms2.scans[linked_msms.scan_index]
                                .precursor_information.charge;

        // We need to ensure that the search succeeded. For example if we are
        // trying to find a number on an empty array, the lower bound will be 0,
        // but that doesn't ensure that we found the right element. A similar
        // situation applies if we try to find a needle not contained on the
        // haystack, in which case the search will return the closest lower
        // bound on the search array or last index of the array.
        auto index = idents_msms_key[i].index;
        auto ident = link_table_idents[index];
        if (linked_msms.msms_id != ident.msms_id) {
            // std::cout << "not found" << std::endl;
            // Generate a theoretical_isotope_distribution based on averagine.
            // auto midas = MIDAs(charge_state = charge_state);
        } else {
            // Generate a theoretical_isotope_distribution based on the given
            // sequence.
            // TODO: Include modifications? There is probably more to it than
            // just calling midas with a sequence and a charge state.
            auto sequence = ident_data.spectrum_ids[ident.entity_id].sequence;
            auto [mzs, perc] =
                theoretical_isotopes_peptide(sequence, charge_state, 0.01);
            auto &peak = peaks[linked_msms.entity_id];

            // We use the retention time of the APEX of the matched peak,
            // not the msms event.
            double peak_rt = peak.local_max_rt;
            double peak_rt_sigma = peak.raw_roi_sigma_mz;
            auto maybe_feature =
                build_feature(peaks_in_use, peaks, peaks_rt_key, mzs, perc,
                              peak_rt_sigma, peak_rt, discrepancy_threshold);
            if (maybe_feature) {
                auto feature = maybe_feature.value();
                features.push_back(feature);
                // Remove used peaks on the feature from the pool.
                for (const auto &peak_id : feature.peak_ids) {
                    peaks_in_use[peak_id] = true;
                }
                // DEBUG: ...
                // std::cout << " rt: " << maybe_feature.value().rt
                //<< " monoisotopic_mz: "
                //<< maybe_feature.value().monoisotopic_mz
                //<< " monoisotopic_height: "
                //<< maybe_feature.value().monoisotopic_height
                //<< " average_mz: " << maybe_feature.value().average_mz
                //<< " total_height: "
                //<< maybe_feature.value().total_height
                //<< " num_isotopes: "
                //<< maybe_feature.value().peak_ids.size() << std::endl;
                // break; // DEBUG: <---
            }
        }
    }
    // TODO: Sort features by height and assign feature ids.
    return features;
}
