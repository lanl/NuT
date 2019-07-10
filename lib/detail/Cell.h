/**
 * Header file for cell class definition
 *
 * @author: Hans R. Hammer
 */

#pragma once

#include "Mesh_Element.h"

namespace nut {

/**
 * @class to describe a geometry cell
 */
struct Cell : public Mesh_Element {

  explicit Cell(index_t id = index_t{}) : Mesh_Element(id) {}
};

} // namespace nut

namespace std {
template<>
struct hash<nut::Cell> {
  using argument_type = nut::Cell;
  using result_type = std::size_t;
  result_type operator()(argument_type const & m) const noexcept {
    return std::hash<nut::index_t>()(m.as_id());
  } // operator()
};
} // namespace std
