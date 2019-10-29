#include "protein_inference/protein_inference.hpp"

#include <algorithm>
#include <map>

// DEBUG:...
#include <iostream>

ProteinInference::Graph ProteinInference::create_graph(
    const IdentData::IdentData &ident_data) {
    std::map<std::string, uint64_t> protein_map;
    std::map<std::string, uint64_t> psm_map;
    ProteinInference::Graph graph;
    // Initialize nodes.
    for (size_t i = 0; i < ident_data.protein_hypotheses.size(); ++i) {
        const auto &protein_hypothesis = ident_data.protein_hypotheses[i];

        // Check if protein already in the graph.
        if (protein_map.count(protein_hypothesis.db_sequence_id) == 0) {
            ProteinInference::Node node = {};
            node.type = ProteinInference::PROTEIN;
            node.id = protein_hypothesis.db_sequence_id;
            node.num = 0;
            node.nodes = {};
            protein_map[protein_hypothesis.db_sequence_id] =
                graph.protein_nodes.size();
            graph.protein_nodes.push_back(node);
        }
        for (const auto &psm : protein_hypothesis.spectrum_ids) {
            if (psm_map.count(psm) == 0) {
                ProteinInference::Node node = {};
                node.type = ProteinInference::PSM;
                node.id = psm;
                node.num = 0;
                node.nodes = {};
                psm_map[psm] = graph.psm_nodes.size();
                graph.psm_nodes.push_back(node);
            }
        }
    }
    // Create graph by assigning pointer to references
    for (size_t i = 0; i < ident_data.protein_hypotheses.size(); ++i) {
        const auto &protein_hypothesis = ident_data.protein_hypotheses[i];

        // Check if protein already in the graph.
        auto protein_index = protein_map.at(protein_hypothesis.db_sequence_id);

        for (const auto &psm : protein_hypothesis.spectrum_ids) {
            auto psm_index = psm_map.at(psm);

            // Update nodes.
            Node *psm_ptr = &graph.psm_nodes[psm_index];
            graph.protein_nodes[protein_index].nodes.push_back(psm_ptr);
            ++graph.protein_nodes[protein_index].num;

            Node *protein_ptr = &graph.protein_nodes[protein_index];
            graph.psm_nodes[psm_index].nodes.push_back(protein_ptr);
            ++graph.psm_nodes[psm_index].num;
        }
    }
    // std::cout << "PROTEINS" << std::endl;
    // for (const auto &node : graph.protein_nodes) {
    // std::cout << "======" << std::endl;
    // std::cout << node.id << " ";
    // std::cout << node.type << " ";
    // std::cout << node.num << " ";
    // std::cout << node.nodes.size() << " ";
    // std::cout << std::endl;
    // std::cout << "------" << std::endl;
    // for (size_t i = 0; i < node.nodes.size(); ++i) {
    // std::cout << node.nodes[i]->id << std::endl;
    //}
    // std::cout << std::endl;
    //}
    // std::cout << "PSM" << std::endl;
    // for (const auto &node : graph.psm_nodes) {
    // std::cout << "======" << std::endl;
    // std::cout << node.id << " ";
    // std::cout << node.type << " ";
    // std::cout << node.num << " ";
    // std::cout << node.nodes.size() << " ";
    // std::cout << std::endl;
    // std::cout << "------" << std::endl;
    // for (size_t i = 0; i < node.nodes.size(); ++i) {
    // std::cout << node.nodes[i]->id << std::endl;
    //}
    // std::cout << std::endl;
    //}
    return graph;
}

void ProteinInference::razor(ProteinInference::Graph &graph) {
    // Sort the protein nodes in descending number of PSM contained in it. We
    // can't modify the graph directly, as this would change the memory pointed
    // by the Node::nodes (Node *).
    //
    // We create a new vector of Node * that we
    // are able to sort while preserving the original order on the graph.
    std::vector<Node *> protein_nodes(graph.protein_nodes.size());
    for (size_t i = 0; i < graph.protein_nodes.size(); ++i) {
        protein_nodes[i] = &graph.protein_nodes[i];
    }
    auto sort_protein_nodes = [](auto a, auto b) { return a->num > b->num; };
    std::stable_sort(protein_nodes.begin(), protein_nodes.end(),
                     sort_protein_nodes);

    for (size_t i = 0; i < protein_nodes.size(); ++i) {
        auto &reference_protein_node = protein_nodes[i];
        for (auto &reference_protein_psm_ptr : reference_protein_node->nodes) {
            if (reference_protein_psm_ptr == nullptr) {
                continue;
            }

            // Visit the other proteins linked to this PSM to remove it's
            // reference.
            for (auto &protein_node_ptr : reference_protein_psm_ptr->nodes) {
                if (protein_node_ptr->id == reference_protein_node->id) {
                    continue;
                }
                for (auto &psm_node_ptr : protein_node_ptr->nodes) {
                    if (psm_node_ptr == nullptr) {
                        continue;
                    }
                    if (psm_node_ptr->id == reference_protein_psm_ptr->id) {
                        psm_node_ptr = nullptr;
                        --protein_node_ptr->num;
                        break;
                    }
                }
            }
        }
        // Sort reminder protein nodes.
        std::stable_sort(protein_nodes.begin() + i, protein_nodes.end(),
                         sort_protein_nodes);
    }
    for (const auto &node : protein_nodes) {
        std::cout << node->id << " " << node->type << " " << node->num
                  << std::endl;
    }
}
