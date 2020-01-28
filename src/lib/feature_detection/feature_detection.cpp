#include "MIDAs/MIDAs.h"

#include "feature_detection/feature_detection.hpp"

FeatureDetection::TheoreticalIsotopes normalize_isotopic_distribution(
    std::vector<Isotopic_Distribution> &isotopes, int8_t charge_state,
    double min_perc) {
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

FeatureDetection::TheoreticalIsotopes
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
    return normalize_isotopic_distribution(isotopes, charge_state, min_perc);
}

// FIXME: This is not entirely accurate. The algorithm is affected by the way
// we do the rounding. For example, MS-Isotope seems to round to the closest
// element instead of the lowest, and then adjusts the number of Hydrogens
// up/down, instead of only up.
struct Element {
    std::string name;
    double proportion;
    double mw;
};
static std::vector<Element> averagine = {
    {"C", 4.9384 / 111.1254, 12.011}, {"H", 7.7583 / 111.1254, 1.008},
    {"N", 1.3577 / 111.1254, 14.007}, {"O", 1.4773 / 111.1254, 15.999},
    {"S", 0.0417 / 111.1254, 32.066},
};
FeatureDetection::TheoreticalIsotopes theoretical_isotopes_formula(
    const std::vector<Element> &average_molecular_composition, double mz,
    int8_t charge_state, double min_perc) {
    // Identify the number of whole atoms that we need. Since there is
    // rounding, we need to adjust by adding hydrogen atoms to be
    // at the same mass. Since we are always rounding down, we will only have
    // mass deficit, not surplus.
    std::vector<Element> num_atoms;
    int32_t hydrogen_index = -1;
    double cumulative_mw = 0.0;
    for (size_t i = 0; i < average_molecular_composition.size(); ++i) {
        const auto &atom = average_molecular_composition[i];
        // Round the number of atoms for this element.
        uint64_t rounded_atoms = atom.proportion * mz * charge_state;
        if (rounded_atoms == 0) {
            continue;
        }
        num_atoms.push_back(
            {atom.name, static_cast<double>(rounded_atoms), atom.mw});
        // Calculate the total molecular weight after the rounding.
        // TODO: Should this lookup be done on a hash table instead of carrying
        // it around in the Element object?
        cumulative_mw += rounded_atoms * atom.mw;
        // Find the index for Hydrogen if we have it.
        if (atom.name == "H") {
            hydrogen_index = i;
        }
    }
    // Calculate the difference between the molecular weight after rounding and
    // the expected mass. The number is calculated in whole units, as this is
    // the number of Hydrogens we are going to add to our table.
    uint64_t extra_hydrogens = mz * charge_state - cumulative_mw;
    if (hydrogen_index == -1) {
        num_atoms.push_back({"H", 0.0, 1.008});
        hydrogen_index = num_atoms.size() - 1;
    }
    num_atoms[hydrogen_index].proportion += extra_hydrogens;

    // Build the chemical formula for the given mz and charge_state.
    std::string formula = "";
    for (const auto &atom : num_atoms) {
        uint64_t num_atoms = atom.proportion;
        formula += atom.name + std::to_string(num_atoms);
    }
    auto midas = MIDAs(charge_state = charge_state);
    auto isotopes = midas.Coarse_Grained_Isotopic_Distribution();
    return normalize_isotopic_distribution(isotopes, charge_state, min_perc);
}

std::optional<FeatureDetection::Feature> FeatureDetection::build_feature(
    const std::vector<bool> &peaks_in_use,
    const std::vector<Centroid::Peak> &peaks,
    const std::vector<Search::KeySort<double>> &peaks_rt_key,
    const TheoreticalIsotopes &theoretical_isotopes, double tolerance_rt,
    double retention_time, double discrepancy_threshold, int8_t charge_state) {
    const auto &mzs = theoretical_isotopes.mzs;
    const auto &perc = theoretical_isotopes.percs;
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
        if ((peak.fitted_mz + peak.fitted_sigma_mz * 2) < mzs[0] ||
            (peak.fitted_mz - peak.fitted_sigma_mz * 2) > mzs[mzs.size() - 1]) {
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
                  return (p1.fitted_mz < p2.fitted_mz);
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
            if ((peak.fitted_mz + peak.fitted_sigma_mz) > theoretical_mz &&
                (peak.fitted_mz - peak.fitted_sigma_mz) < theoretical_mz) {
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
                double normalized_height =
                    candidate->fitted_height / ref_candidate->fitted_height;
                double a = candidate->fitted_mz - theoretical_mz;
                double b = candidate->fitted_rt - retention_time;
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
    feature.msms_id = -1;
    // FIXME: Currently assuming that the minimum detected isotope is the
    // monoisotopic peak, but THIS MIGHT NOT BE THE CASE. For simplicity and to
    // keep the flow going I'll leave this for now, but must go back and FIX
    // it.
    feature.monoisotopic_mz = selected_candidates[0]->fitted_mz;
    feature.monoisotopic_height = selected_candidates[0]->fitted_height;
    // Find the weighted average mz and the average height of the selected
    // isotopes.
    feature.average_mz = 0.0;
    feature.total_height = 0.0;
    feature.average_rt = 0.0;
    feature.average_rt_delta = 0.0;
    feature.average_rt_sigma = 0.0;
    feature.average_mz_sigma = 0.0;
    for (size_t i = 0; i < selected_candidates.size(); ++i) {
        auto candidate = selected_candidates[i];
        feature.total_height += candidate->fitted_height;
        feature.average_mz += candidate->fitted_height * candidate->fitted_mz;
        feature.average_rt += candidate->fitted_rt;
        feature.average_rt_delta += candidate->rt_delta;
        feature.average_rt_sigma += candidate->fitted_sigma_rt;
        feature.average_mz_sigma += candidate->fitted_sigma_mz;
        feature.peak_ids.push_back(candidate->id);
    }
    if (feature.total_height == 0) {
        return std::nullopt;
    }
    feature.average_mz /= feature.total_height;
    feature.average_rt /= selected_candidates.size();
    feature.average_rt_delta /= selected_candidates.size();
    feature.average_rt_sigma /= selected_candidates.size();
    feature.average_mz_sigma /= selected_candidates.size();
    feature.charge_state = charge_state;
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
    // 0.- Copy and sort the necessary vectors for the use of binary search
    //     (link_table_idents is sorted by msms, as that is the key being
    //     used for searching). The peaks array and link_table_msms array
    //     should have already been sorted by peak_id, which is what we
    //     want, as we are going to start matching peaks from highest
    //     intensity to lowest.
    // 1.- For each linked peak on the link_table_msms, find it's associated
    //     entry on the link_table_idents.
    // 2.- If the entry is found, use it to generate a theoretical isotopic
    //     distribution, otherwise, averagine will be generated.
    // 3.- We try to find the proposed peaks from the theoretical
    //     distribution on the peaks array (Maximum likelihood).
    // 4.- The peaks are marked as non available for future use. This is
    //     a greedy algorithm.

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
        peaks_rt_key[i] = {i, peaks[i].fitted_rt};
    }
    {
        auto sorting_key_func = [](const Search::KeySort<double> &p1,
                                   const Search::KeySort<double> &p2) {
            return (p1.sorting_key < p2.sorting_key);
        };
        std::sort(peaks_rt_key.begin(), peaks_rt_key.end(), sorting_key_func);
    }
    // We use this variable to keep track of the peaks we have already linked.
    auto peaks_in_use = std::vector<bool>(peaks.size());

    // TODO: We should probably prioritize the MSMS events that HAVE an
    // identification, instead of being intensity based only.
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
        const auto &peak = peaks[linked_msms.entity_id];
        TheoreticalIsotopes theoretical_isotopes = {};
        if (linked_msms.msms_id != ident.msms_id) {
            // Generate a theoretical_isotope_distribution based on averagine.
            // FIXME: This needs to be checked. For now we are ignoring the
            // msms events without identifications to get initial measurements
            // with our data.
            // theoretical_isotopes = theoretical_isotopes_formula(
            // averagine, peak.fitted_mz, charge_state, 0.01);
        } else {
            // Generate a theoretical_isotope_distribution based on the
            // given sequence.
            // FIXME: This algorithm WILL be affected by modifications. If we
            // don't consider them, the m/z results will not match the peaks in
            // our peak list.
            auto sequence = ident_data.spectrum_ids[ident.entity_id].sequence;
            theoretical_isotopes =
                FeatureDetection::theoretical_isotopes_peptide(
                    sequence, charge_state, 0.01);
        }
        if (theoretical_isotopes.mzs.empty() ||
            theoretical_isotopes.percs.empty()) {
            continue;
        }

        // We use the retention time of the APEX of the matched peak,
        // not the msms event.
        double peak_rt = peak.fitted_rt;
        double peak_rt_sigma = peak.fitted_sigma_rt;
        auto maybe_feature = build_feature(
            peaks_in_use, peaks, peaks_rt_key, theoretical_isotopes,
            peak_rt_sigma * 2, peak_rt, discrepancy_threshold, charge_state);
        if (maybe_feature) {
            auto feature = maybe_feature.value();
            feature.msms_id = linked_msms.msms_id;
            features.push_back(feature);
            // Remove used peaks on the feature from the pool.
            for (const auto &peak_id : feature.peak_ids) {
                peaks_in_use[peak_id] = true;
            }
        }
    }

    // Sort features by height and assign feature ids.
    std::sort(features.begin(), features.end(), [](auto &a, auto &b) {
        return (a.total_height > b.total_height);
    });
    for (size_t i = 0; i < features.size(); ++i) {
        features[i].id = i;
    }

    return features;
}

#include <iostream>

std::vector<std::vector<uint64_t>> find_all_paths(
    FeatureDetection::CandidateGraph &graph, uint64_t root_node) {
    std::vector<std::vector<uint64_t>> paths;
    std::vector<std::vector<uint64_t>> stack;
    stack.push_back({root_node});
    while (!stack.empty()) {
        auto curr_path = stack.back();
        stack.pop_back();
        auto last = curr_path[curr_path.size() - 1];
        auto &root_node = graph[last];
        if (root_node.visited) {
            continue;
        }

        bool path_finished = true;
        for (const auto &node : root_node.nodes) {
            if (graph[node].visited) {
                continue;
            }
            path_finished = false;
            auto new_path = curr_path;
            new_path.push_back(node);
            stack.push_back(new_path);
        }
        if (path_finished) {
            paths.push_back(curr_path);
        }
    }

    return paths;
}

struct RollingCosineResults {
    double best_dot;
    size_t best_shift;
    size_t pad;
};
// NOTE: The order matters. A should be the path we are exploring, and B the
// reference theoretical path.
RollingCosineResults rolling_cosine_sim(std::vector<double> &A,
                                        std::vector<double> &B) {
    // We need at least 2 points to form a feature.
    if (A.size() < 2 || B.size() < 2) {
        return {0.0, 0, 0};
    }
    // Pre-calculate the norm of A and B.
    double norm_a = 0.0;
    for (size_t i = 0; i < A.size(); ++i) {
        norm_a += A[i] * A[i];
    }
    norm_a = std::sqrt(norm_a);
    double norm_b = 0.0;
    double max_b = 0.0;
    size_t max_b_index = 0;
    for (size_t i = 0; i < B.size(); ++i) {
        norm_b += B[i] * B[i];
        if (B[i] > max_b) {
            max_b = B[i];
            max_b_index = i;
        }
    }
    norm_b = std::sqrt(norm_b);
    double denom = norm_a * norm_b;
    // Create a left padded version of A. This is used to shift B over A for
    // the calculation of cosine similarity. This shift ensures that the
    // maximum in sequence B is tested in all positions of A.
    size_t pad = max_b_index;
    std::vector<double> C = std::vector<double>(pad + A.size(), 0.0);
    for (size_t i = 0; i < A.size(); ++i) {
        C[i + pad] = A[i];
    }
    // Calculate the dot product for all shifts of A and keep the maximum.
    double best_dot = 0.0;
    size_t best_shift = 0;  // `shift` tell us where in A we have the max of B.
    for (size_t k = 0; k < C.size() - 1; ++k) {
        double dot = 0.0;
        for (size_t i = 0; i < std::max(A.size(), B.size()); ++i) {
            double a = 0;
            double b = 0;
            if (i + k < C.size()) {
                a = C[i + k];
            }
            if (i < B.size()) {
                b = B[i];
            }
            dot += a * b;
        }
        dot = dot / denom;
        if (dot > best_dot) {
            best_dot = dot;
            best_shift = k;
        }
    }
    return {best_dot, best_shift, pad};
}

// Trim path A based on the results of rolling cosine similarity.
std::vector<size_t> trim_path(const std::vector<size_t> &path,
                              std::vector<double> &A, std::vector<double> &B,
                              RollingCosineResults &sim) {
    // std::vector<size_t> path;
    // std::vector<size_t> trimmed_path;
    // for (size_t i = 0; i < B.size(); ++i) {
    //}
    //
    // if shift < pad (we are on the left side)...
    //
    // x x x x x 0 1 2 3 7 8 x x x x x
    // x x x x A B C D E F x x x x x x
    //
    // ...we can trim right if B.size() + shift  > pad + A.size():
    //
    // x x x x x 0 1 2 3 7 x x x x x x
    // x x x x A B C D E F x x x x x x
    //
    // if shift > pad (we are on the right side)...
    //
    // x x x x x 0 1 2 3 7 8 x x x x x
    // x x x x x x x x A B C D E F x x
    //
    // ...we can trim left:
    //
    // x x x x x x x x 3 7 8 x x x x x
    // x x x x x x x x A B C D E F x x
    //
    // In addition to this, we need to check if there are any discrepancies
    // larger than a given threshold from the theoretical %.
    // NOTE: min_i and max_i are given by the trimming from the previous ASCII
    // scheme.
    size_t min_i = sim.best_shift - sim.pad;
    // size_t max_i = std::min(B.size() - min_i, A.size());
    // size_t min_i = 0;
    size_t max_i = A.size();
    std::vector<size_t> trimmed_path;
    std::cout << "LETS GOOOOOO:" << std::endl;
    std::cout << "best_dot: " << sim.best_dot << std::endl;
    std::cout << "best_shift: " << sim.best_shift << std::endl;
    std::cout << "pad: " << sim.pad << std::endl;
    // std::cout << "charge_state: " << sim.charge_state << std::endl;
    for (size_t i = min_i; i < max_i; ++i) {  // TODO: max_i?
        // Normalize A to the maximum of B (shift).
        double norm_a = A[i] / A[sim.best_shift];
        // Check the difference between the theoretical and normalized
        // intensity.
        double abs_diff = std::abs(norm_a - B[sim.best_shift - sim.pad + i]);
        trimmed_path.push_back(path[i]);
    }
    // return path;
    return trimmed_path;
}

#include <map>
std::map<double, std::vector<double>> averagine_table = {
    {300.17921, {100.00, 16.53, 2.10, 0.20, 0.01}},
    {400.0, {100.00, 22.68, 3.48, 0.40, 0.04}},
    {500.29530, {100.00, 27.70, 5.12, 0.71, 0.08, 0.01}},
    {700.43520, {100.00, 39.22, 9.33, 1.65, 0.24, 0.03}},
    {800.0, {100.00, 45.17, 12.21, 2.45, 0.40, 0.05, 0.01}},
    {999.71361, {100.00, 55.58, 17.81, 4.18, 0.79, 0.13, 0.02}},
    {1999.93307,
     {89.11, 100.00, 64.45, 30.32, 11.43, 3.62, 1.00, 0.24, 0.05, 0.01}},
};
void FeatureDetection::find_candidates(
    const std::vector<Centroid::Peak> &peaks,
    const std::vector<uint8_t> &charge_states) {
    double carbon_diff = 1.0033;  // NOTE: Maxquant uses 1.00286864

    // Sort peaks by mz.
    auto sorted_peaks = std::vector<Search::KeySort<double>>(peaks.size());
    for (size_t i = 0; i < peaks.size(); ++i) {
        sorted_peaks[i] = {i, peaks[i].fitted_mz};
    }
    std::stable_sort(
        sorted_peaks.begin(), sorted_peaks.end(),
        [](auto &p1, auto &p2) { return (p1.sorting_key < p2.sorting_key); });

    // Initialize graph.
    std::vector<FeatureDetection::CandidateGraph> charge_state_graphs(
        charge_states.size());
    for (size_t k = 0; k < charge_states.size(); ++k) {
        charge_state_graphs[k].resize(sorted_peaks.size());
    }
    for (size_t i = 0; i < sorted_peaks.size(); ++i) {
        auto &ref_peak = peaks[sorted_peaks[i].index];
        // DEBUG: Change to the commented version when ready.
        double tol_mz = 0.1;
        double tol_rt = 0.1;
        // double tol_mz = ref_peak.fitted_sigma_mz;
        // double tol_rt = ref_peak.fitted_sigma_rt;
        double min_rt = ref_peak.fitted_rt - tol_rt;
        double max_rt = ref_peak.fitted_rt + tol_rt;
        for (size_t k = 0; k < charge_states.size(); ++k) {
            charge_state_graphs[k][i].id = ref_peak.id;
            auto charge_state = charge_states[k];
            if (charge_state == 0) {
                continue;
            }
            double mz_diff = carbon_diff / charge_state;
            double min_mz = (ref_peak.fitted_mz + mz_diff) - tol_mz;
            double max_mz = (ref_peak.fitted_mz + mz_diff) + tol_mz;

            // Find peaks within tolerance range and add them to the graph.
            for (size_t j = i + 1; j < sorted_peaks.size(); ++j) {
                auto &peak = peaks[sorted_peaks[j].index];
                if (peak.fitted_mz > max_mz) {
                    break;
                }
                if (peak.fitted_mz > min_mz && peak.fitted_rt > min_rt &&
                    peak.fitted_rt < max_rt) {
                    charge_state_graphs[k][i].nodes.push_back(j);
                }
            }
        }
    }

    // Visit nodes to find most likely features.
    size_t i = 0;
    while (i < sorted_peaks.size()) {
        double best_dot = 0.0;
        std::vector<uint64_t> best_path;
        uint8_t best_charge_state = 0;
        for (size_t k = 0; k < charge_states.size(); ++k) {
            int64_t charge_state = charge_states[k];
            auto paths = find_all_paths(charge_state_graphs[k], i);
            double ref_mz =
                peaks[sorted_peaks[i].index].fitted_mz * charge_states[k];
            auto ref_key_it = averagine_table.lower_bound(ref_mz);
            if (ref_key_it == averagine_table.end()) {
                continue;
            }
            // Check if the next element is closer to the ref_mz.
            auto ref_key_it_next = ref_key_it;
            ref_key_it_next++;
            if (ref_key_it_next != averagine_table.end() &&
                ref_mz - ref_key_it->first > ref_key_it_next->first - ref_mz) {
                ref_key_it++;
            }
            std::vector<double> &ref_heights = ref_key_it->second;
            for (const auto &path : paths) {
                if (path.size() < 2) {
                    continue;
                }
                std::vector<double> path_heights;
                for (const auto &x : path) {
                    path_heights.push_back(
                        peaks[sorted_peaks[x].index].fitted_height);
                }
                auto sim = rolling_cosine_sim(path_heights, ref_heights);
                if (sim.best_dot > best_dot) {
                    best_dot = sim.best_dot;
                    best_charge_state = charge_state;
                    best_path = trim_path(path, path_heights, ref_heights, sim);
                }
            }
        }
        if (best_path.empty()) {
            ++i;
            continue;
        }
        // We found a feature, but the initial reference is not included.
        if (best_path[0] == i) {
            ++i;
        }
        std::cout << "BEST DOT: " << best_dot << " for PATH: ";
        for (const auto &x : best_path) {
            std::cout << peaks[sorted_peaks[x].index].id << ' ';
            // TODO: Build feature.
            // Mark peaks as used.
            for (size_t k = 0; k < charge_states.size(); ++k) {
                charge_state_graphs[k][x].visited = true;
            }
        }
        std::cout << " Z: " << (int64_t)best_charge_state << std::endl;
        std::cout << std::endl;

        // break;
    }
}
