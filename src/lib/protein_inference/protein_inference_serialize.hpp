#ifndef PROTEININFERENCE_PROTEININFERENCE_SERIALIZE_HPP
#define PROTEININFERENCE_PROTEININFERENCE_SERIALIZE_HPP

#include <iostream>

#include "protein_inference/protein_inference.hpp"

// This namespace groups the functions used to serialize ProteinInference data
// structures into a binary stream.
namespace ProteinInference::Serialize {

bool read_inferred_protein(std::istream &stream,
                           ProteinInference::InferredProtein *protein);
bool write_inferred_protein(std::ostream &stream,
                            const ProteinInference::InferredProtein &protein);

bool read_inferred_proteins(
    std::istream &stream,
    std::vector<ProteinInference::InferredProtein> *proteins);
bool write_inferred_proteins(std::ostream &stream,
                             const std::vector<InferredProtein> &proteins);

}  // namespace ProteinInference::Serialize

#endif /* PROTEININFERENCE_PROTEININFERENCE_SERIALIZE_HPP */
