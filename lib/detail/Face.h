/**
 * Header file for the surface class
 *
 * @author: Hans R. Hammer
 *
 */

#pragma once

#include "Mesh_Element.h"
#include <functional> // hash

namespace nut {

struct Face : public Mesh_Element {

  explicit Face(index_t id) : Mesh_Element(id) {}
};

struct Cartesian_Face : public Mesh_Element {
  Cartesian_Face(index_t index) : Mesh_Element(index) {}

  index_t get_index() const {
    return (as_id());
  }
};

struct Spherical_1D_Face : public Mesh_Element {
  Spherical_1D_Face(index_t index) : Mesh_Element(index) {}

  index_t get_index() const {
    return as_id();
  }
};

} // namespace nut

namespace std {
template<>
struct hash<nut::Spherical_1D_Face> {
  using argument_type = nut::Spherical_1D_Face;
  using result_type = std::size_t;
  result_type operator()(argument_type const & m) const noexcept {
    return std::hash<nut::index_t>()(m.as_id());
  } // operator()
};

template<>
struct hash<nut::Cartesian_Face> {
  using argument_type = nut::Cartesian_Face;
  using result_type = std::size_t;
  result_type operator()(argument_type const & m) const noexcept {
    return std::hash<nut::index_t>()(m.as_id());
  } // operator()
};
} // namespace std

// End of file
