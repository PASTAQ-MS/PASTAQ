#include "protein_inference/protein_inference_serialize.hpp"
#include "utils/serialization.hpp"

bool ProteinInference::Serialize::read_inferred_protein(
    std::istream &stream, ProteinInference::InferredProtein *protein) {
    Serialization::read_string(stream, &protein->protein_id);
    Serialization::read_string(stream, &protein->psm_id);
    return stream.good();
}

bool ProteinInference::Serialize::write_inferred_protein(
    std::ostream &stream, const ProteinInference::InferredProtein &protein) {
    Serialization::write_string(stream, protein.protein_id);
    Serialization::write_string(stream, protein.psm_id);
    return stream.good();
}

bool ProteinInference::Serialize::read_inferred_proteins(
    std::istream &stream,
    std::vector<ProteinInference::InferredProtein> *proteins) {
    Serialization::read_vector<ProteinInference::InferredProtein>(
        stream, proteins, read_inferred_protein);
    return stream.good();
}

bool ProteinInference::Serialize::write_inferred_proteins(
    std::ostream &stream,
    const std::vector<ProteinInference::InferredProtein> &proteins) {
    Serialization::write_vector<ProteinInference::InferredProtein>(
        stream, proteins, write_inferred_protein);
    return stream.good();
}
