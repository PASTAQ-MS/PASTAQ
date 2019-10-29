#ifndef PROTEININFERENCE_PROTEININFERENCE_HPP
#define PROTEININFERENCE_PROTEININFERENCE_HPP

#include <optional>

#include "raw_data/raw_data.hpp"

namespace ProteinInference {
// This is a node for the bipartite graph representing the relationships
// between protein hypotheses and peptide spectrum matches (PSM).
enum NodeType : uint8_t { UNKNOWN = 0, PROTEIN = 1, PSM = 2 };
// TODO: Explain why this structure is the way it is.
struct Node {
    NodeType type;
    uint64_t num;
    std::string id;
    std::vector<std::optional<uint64_t>> nodes;
};

struct Graph {
    std::vector<Node> protein_nodes;
    std::vector<Node> psm_nodes;
};

// Initializes the initial graph. The returning vector corresponds to the
// protein nodes, which contain links to their respective PSMs.
Graph create_graph(const IdentData::IdentData &ident_data);

// Performs Occam's razor protein inference, where we select the minimum number
// of proteins that explain the observed PSM in the graph. Note that the graph
// will be modified in place.
void razor(Graph &graph);

}  // namespace ProteinInference

#endif /* PROTEININFERENCE_PROTEININFERENCE_HPP */
