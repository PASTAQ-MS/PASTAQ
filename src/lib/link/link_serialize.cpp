#include "link/link_serialize.hpp"
#include "utils/serialization.hpp"

bool Link::Serialize::read_linked_msms(std::istream &stream,
                                       Link::LinkedMsms *link) {
    Serialization::read_uint64(stream, &link->entity_id);
    Serialization::read_uint64(stream, &link->msms_id);
    Serialization::read_uint64(stream, &link->scan_index);
    Serialization::read_double(stream, &link->distance);
    return stream.good();
}

bool Link::Serialize::write_linked_msms(std::ostream &stream,
                                        const Link::LinkedMsms &link) {
    Serialization::write_uint64(stream, link.entity_id);
    Serialization::write_uint64(stream, link.msms_id);
    Serialization::write_uint64(stream, link.scan_index);
    Serialization::write_double(stream, link.distance);
    return stream.good();
}

bool Link::Serialize::read_linked_msms_table(
    std::istream &stream, std::vector<Link::LinkedMsms> *links) {
    uint64_t num_rows = 0;
    Serialization::read_uint64(stream, &num_rows);
    *links = std::vector<Link::LinkedMsms>(num_rows);
    for (size_t i = 0; i < num_rows; ++i) {
        Link::Serialize::read_linked_msms(stream, &(*links)[i]);
    }
    return stream.good();
}

bool Link::Serialize::write_linked_msms_table(
    std::ostream &stream, const std::vector<Link::LinkedMsms> &links) {
    uint64_t num_rows = links.size();
    Serialization::write_uint64(stream, num_rows);
    for (size_t i = 0; i < num_rows; ++i) {
        Link::Serialize::write_linked_msms(stream, links[i]);
    }
    return stream.good();
}

