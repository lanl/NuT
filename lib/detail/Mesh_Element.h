// mesh_element.h
// Jan 08, 2019
// (c) Copyright 2019 LANSLLC, all rights reserved

#pragma once

#include "types.hh"
#include <functional>  // hash

namespace nut {

struct Mesh_Element {
public:
  explicit Mesh_Element(index_t const id) : id_(id) {}

  bool operator==(Mesh_Element const & e) const { return e.id_ == id_; }

  bool operator!=(Mesh_Element const & e) const { return e.id_ != id_; }

  bool operator<(Mesh_Element const & other) const
  {
    return this->id_ < other.id_;
  }

  /* Had to rename id() method. */
  index_t as_id() const { return id_; }

private:
  index_t id() const { return id_; }

  index_t id_;

};  // mesh_element

// inline bool
// operator<(Mesh_Element const & m1, Mesh_Element const & m2) {
//   return m1.id() < m2.id();
// }

struct Mesh_Element_Hash {
  std::size_t operator()(Mesh_Element const & m) const noexcept
  {
    return std::hash<index_t>()(m.as_id());
  }
};

}  // namespace nut

namespace std {
template <>
struct hash<nut::Mesh_Element> {
  using argument_type = nut::Mesh_Element;
  using result_type = std::size_t;
  result_type operator()(argument_type const & m) const noexcept
  {
    return std::hash<nut::index_t>()(m.as_id());
  }  // operator()
};
}  // namespace std

// End of file
