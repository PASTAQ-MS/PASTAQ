#ifndef PROTEININFERENCE_PROTEININFERENCE_HPP
#define PROTEININFERENCE_PROTEININFERENCE_HPP

#include <optional>

#include "raw_data/raw_data.hpp"

namespace ProteinInference {
enum NodeType : uint8_t { UNKNOWN = 0, PROTEIN = 1, PSM = 2 };

// This is a node for the bipartite graph representing the relationships
// between protein hypotheses and peptide spectrum matches (PSM).
//
// Each Node contains information regarding it's type, the number of elements
// it contains and a unique identifier. The identifier must be unique within
// each NodeType. For example, we can use the same identifier in a PROTEIN node
// and a PSM node.
//
// The graph is represented using adjacency lists, as the protein inference is
// a highly sparse graph. The adjacency list on each node consists on the index
// of a node for the opposite node type, i.e., a PROTEIN node's adjacency list
// points to the indices in the PSM vector for the graph.
//
// For performance reasons, we allow the adjacency list nodes to be null. This
// allow us to cut a link between PROTEIN/PSM without having to regenerate
// the adjacency list or copy memory around.
struct Node {
    NodeType type;
    uint64_t num;
    std::string id;
    std::vector<std::optional<uint64_t>> nodes;
};

// This structure owns the memory for the different nodes in the protein
// inference problem. The nodes of each type are stored in a different vector.
// Note that due to the nature of the Node data structure, we might have
// unlinked nodes in the graph after performing protein inference. This is by
// design, and if need be, in the future we can use a reduction step to remove
// unused nodes and connections.
struct Graph {
    std::vector<Node> protein_nodes;
    std::vector<Node> psm_nodes;
};

// TODO: Is this really the best way of associating PSM with Proteins? Should
// we reference the index of the ident files?
// TODO: Add documentation for this.
struct InferredProtein {
    std::string protein_id;
    std::string psm_id;
};

// Initializes the initial graph. The returning vector corresponds to the
// protein nodes, which contain links to their respective PSMs.
Graph create_graph(const IdentData::IdentData &ident_data);

// Performs Occam's razor protein inference, where we select the minimum number
// of proteins that explain the observed PSM in the graph. Note that the graph
// will be modified in place.
std::vector<InferredProtein> razor(const IdentData::IdentData &ident_data);

}  // namespace ProteinInference

#endif /* PROTEININFERENCE_PROTEININFERENCE_HPP */
