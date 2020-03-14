#include "protein_inference/protein_inference.hpp"

#include <algorithm>
#include <map>

ProteinInference::Graph ProteinInference::create_graph(
    const IdentData::IdentData &ident_data) {
    std::map<std::string, uint64_t> protein_map;
    std::map<std::string, uint64_t> psm_map;
    ProteinInference::Graph graph;
    // FIXME: Rework this to use the new IdentData data structures.
    // Initialize nodes.
    // for (size_t i = 0; i < ident_data.protein_hypotheses.size(); ++i) {
    //     const auto &protein_hypothesis = ident_data.protein_hypotheses[i];
    //     if (protein_map.count(protein_hypothesis.db_sequence_id) == 0) {
    //         ProteinInference::Node node = {};
    //         node.type = ProteinInference::PROTEIN;
    //         node.id = protein_hypothesis.db_sequence_id;
    //         node.num = 0;
    //         node.nodes = {};
    //         protein_map[protein_hypothesis.db_sequence_id] =
    //             graph.protein_nodes.size();
    //         graph.protein_nodes.push_back(node);
    //     }
    //     for (const auto &psm : protein_hypothesis.spectrum_ids) {
    //         if (psm_map.count(psm) == 0) {
    //             ProteinInference::Node node = {};
    //             node.type = ProteinInference::PSM;
    //             node.id = psm;
    //             node.num = 0;
    //             node.nodes = {};
    //             psm_map[psm] = graph.psm_nodes.size();
    //             graph.psm_nodes.push_back(node);
    //         }
    //     }
    // }
    // // Create graph by assigning edges to the adjacency list of each node.
    // for (size_t i = 0; i < ident_data.protein_hypotheses.size(); ++i) {
    //     const auto &protein_hypothesis = ident_data.protein_hypotheses[i];

    //     // Check if protein already in the graph.
    //     auto protein_index = protein_map.at(protein_hypothesis.db_sequence_id);

    //     for (const auto &psm : protein_hypothesis.spectrum_ids) {
    //         auto psm_index = psm_map.at(psm);

    //         // Update nodes. We need to make sure that the node was not
    //         // included already on the adjacency list.
    //         auto ptr_in_node_list =
    //             [](uint64_t target_node,
    //                std::vector<std::optional<uint64_t>> &node_list) {
    //                 for (const auto &node : node_list) {
    //                     if (node == target_node) {
    //                         return true;
    //                     }
    //                 }
    //                 return false;
    //             };
    //         if (!ptr_in_node_list(psm_index,
    //                               graph.protein_nodes[protein_index].nodes)) {
    //             graph.protein_nodes[protein_index].nodes.push_back(psm_index);
    //             ++graph.protein_nodes[protein_index].num;
    //         }
    //         if (!ptr_in_node_list(protein_index,
    //                               graph.psm_nodes[psm_index].nodes)) {
    //             graph.psm_nodes[psm_index].nodes.push_back(protein_index);
    //             ++graph.psm_nodes[psm_index].num;
    //         }
    //     }
    // }
    return graph;
}

std::vector<ProteinInference::InferredProtein> ProteinInference::razor(
    const IdentData::IdentData &ident_data) {
    // Initialize the base inference graph.
    auto graph = ProteinInference::create_graph(ident_data);

    // Sort the protein nodes in descending number of PSM contained in it. We
    // can't modify the graph directly, as this would change the indexes
    // referenced in the adjacency lists.
    //
    // We create a new vector of Node pointers that we are able to sort while
    // preserving the original order on the graph.
    std::vector<Node *> protein_nodes(graph.protein_nodes.size());
    for (size_t i = 0; i < graph.protein_nodes.size(); ++i) {
        protein_nodes[i] = &graph.protein_nodes[i];
    }
    auto sort_protein_nodes = [](auto a, auto b) { return a->num > b->num; };
    std::stable_sort(protein_nodes.begin(), protein_nodes.end(),
                     sort_protein_nodes);

    for (size_t i = 0; i < protein_nodes.size(); ++i) {
        // This is the protein at the top of the list, i.e., the protein that
        // is associated with the most peptides.
        auto &ref_protein_ptr = protein_nodes[i];
        for (auto &ref_psm_index : ref_protein_ptr->nodes) {
            if (!ref_psm_index) {
                continue;
            }
            auto &ref_psm = graph.psm_nodes[ref_psm_index.value()];

            // Visit the other proteins linked to this PSM to sever the links
            // from the adjacency lists.
            for (auto &cur_protein_index : ref_psm.nodes) {
                if (!cur_protein_index) {
                    continue;
                }
                auto &cur_protein =
                    graph.protein_nodes[cur_protein_index.value()];
                if (&cur_protein == ref_protein_ptr) {
                    continue;
                }
                for (auto &cur_psm_index : cur_protein.nodes) {
                    if (!cur_psm_index) {
                        continue;
                    }
                    auto &cur_psm = graph.psm_nodes[cur_psm_index.value()];
                    if (&cur_psm == &ref_psm) {
                        --cur_protein.num;
                        --cur_psm.num;
                        cur_psm_index = std::nullopt;
                        cur_protein_index = std::nullopt;
                        break;
                    }
                }
            }
        }
        // Sort reminder protein nodes.
        std::stable_sort(protein_nodes.begin() + i, protein_nodes.end(),
                         sort_protein_nodes);
    }

    // Fill up the result table.
    std::vector<InferredProtein> inferred_proteins;
    for (auto &protein : protein_nodes) {
        for (const auto &psm_index : protein->nodes) {
            if (!psm_index) {
                continue;
            }
            const auto &psm = graph.psm_nodes[psm_index.value()];
            inferred_proteins.push_back({protein->id, psm.id});
        }
    }
    return inferred_proteins;
}
