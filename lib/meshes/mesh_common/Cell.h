/**
 * Header file for cell class definition
 *
 * @author: Hans R. Hammer
 */

#pragma once

#include "mesh_common/Mesh_Element.h"

namespace murmeln_mesh {

/**
 * @class to describe a geometry cell
 */
struct Cell : public Mesh_Element {

  explicit Cell(index_t id = index_t{}) : Mesh_Element(id) {}
};

} // namespace murmeln_mesh

namespace std {
template <> struct hash<murmeln_mesh::Cell> {
  using argument_type = murmeln_mesh::Cell;
  using result_type = std::size_t;
  result_type operator()(argument_type const &m) const noexcept {
    return std::hash<murmeln_mesh::index_t>()(m.as_id());
  } // operator()
};
} // namespace std
