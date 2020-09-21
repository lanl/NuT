// Cell_Face_Descriptor.h
// T. M. Kelley
// Jun 25, 2019
// (c) Copyright 2019 Triad National Security, all rights reserved

#pragma once

#include "mesh_common/boundary.h"

namespace nut_mesh {

/**\brief A data structure that represents the boundaries, if any, on each of
 * six faces. */
template <typename Face_Name, size_t max_faces = 6 /**/>
struct Cell_Face_Descriptor {
public:
  // void set_boundary(Face_Name const & f, boundary::BDY_TYPE const & b);

  void set_boundary(Face_Name const &f, boundary::BDY_TYPE const &b);

  void clear_boundary(Face_Name const &f); // {

  /**\brief Is this cell any kind of boundary? */
  bool is_boundary() const; // {

  /**\brief Is this face any kind of boundary? */
  bool is_boundary(Face_Name const &f) const;

  /**\brief Is the face f set to boundary b? */
  bool is_boundary(Face_Name const &f, boundary::BDY_TYPE const &b) const;

  /**\brief Is the face f set to boundary b? */
  boundary::BDY_TYPE get_boundary(Face_Name const &f) const;

private:
  /* TO DO: this is specific to cubic meshes */
  static constexpr uint32_t clear_face_masks[6] = {
      0b00'000'000'000'000'111'111'111'111'111'000,
      0b00'000'000'000'000'111'111'111'111'000'111,
      0b00'000'000'000'000'111'111'111'000'111'111,
      0b00'000'000'000'000'111'111'000'111'111'111,
      0b00'000'000'000'000'111'000'111'111'111'111,
      0b00'000'000'000'000'000'111'111'111'111'111,
  };

  uint32_t d_{0};
}; // struct Cell_Face_Descriptor

// ---------------- Cell_Face_Descriptor methods ----------------
template <typename Face_Name, size_t mxfcs>
inline void Cell_Face_Descriptor<Face_Name, mxfcs>::set_boundary(
    Face_Name const &f, boundary::BDY_TYPE const &b) {
  if (is_boundary(f, b)) {
    clear_boundary(f);
  }
  // 3 bits per face
  index_t face_idx = (b & 7) << (3 * f);
  d_ |= face_idx;
  return;
}

template <typename Face_Name, size_t mxfcs>
inline void
Cell_Face_Descriptor<Face_Name, mxfcs>::clear_boundary(Face_Name const &f) {
  d_ &= clear_face_masks[f];
  return;
}

template <typename Face_Name, size_t mxfcs>
inline bool Cell_Face_Descriptor<Face_Name, mxfcs>::is_boundary() const {
  return d_ != boundary::NONE;
}

template <typename Face_Name, size_t mxfcs>
inline bool
Cell_Face_Descriptor<Face_Name, mxfcs>::is_boundary(Face_Name const &f) const {
  boundary::BDY_TYPE fbase{get_boundary(f)};
  return fbase != boundary::NONE;
}

template <typename Face_Name, size_t mxfcs>
inline bool Cell_Face_Descriptor<Face_Name, mxfcs>::is_boundary(
    Face_Name const &f, boundary::BDY_TYPE const &b) const {
  boundary::BDY_TYPE fbase{get_boundary(f)};
  return fbase == b;
}

template <typename Face_Name, size_t mxfcs>
inline boundary::BDY_TYPE
Cell_Face_Descriptor<Face_Name, mxfcs>::get_boundary(Face_Name const &f) const {
  index_t const shifter = 3 * f;
  index_t fbase = ((static_cast<index_t>(d_)) >> shifter) & 7;
  return boundary::BDY_TYPE{fbase};
}

} // namespace nut_mesh

// End of file
