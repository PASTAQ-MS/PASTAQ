#ifndef LINK_LINKMSMS_SERIALIZE_HPP
#define LINK_LINKMSMS_SERIALIZE_HPP

#include <iostream>

#include "link/link.hpp"

// This namespace groups the functions used to serialize Link data structures
// into a binary stream.
namespace Link::Serialize {
// TODO(alex): This needs more documentation
bool read_linked_msms(std::istream &stream, Link::LinkedMsms *link);
bool write_linked_msms(std::ostream &stream, const Link::LinkedMsms &link);
bool read_linked_msms_table(std::istream &stream,
                            std::vector<Link::LinkedMsms> *links);
bool write_linked_msms_table(std::ostream &stream,
                             const std::vector<LinkedMsms> &links);
bool read_linked_psm_table(std::istream &stream,
                           std::vector<Link::LinkedPsm> *links);
bool write_linked_psm_table(std::ostream &stream,
                            const std::vector<LinkedPsm> &links);

}  // namespace Link::Serialize

#endif /* LINK_LINKMSMS_SERIALIZE_HPP */
