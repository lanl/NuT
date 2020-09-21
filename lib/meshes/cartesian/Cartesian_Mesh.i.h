// Cartesian_Mesh.i.hh
// May 21, 2019
// (c) Copyright 2019 Triad National Security, all rights reserved

#pragma once

// Cartesian_Mesh inline implementations

//
#ifndef IM_ALLOWED_TO_INCLUDE_CARTESIAN_MESH_I_H
#error "Use Cartesian_Mesh.hh---do not include Cartesian_Mesh.i.hh"
#endif

#include <cmath>
#include <iterator>
#include <string>

namespace nut_mesh {

// ---------------- Cartesian_Mesh methods ----------------
inline Cartesian_Mesh::Vector Cartesian_Mesh::reflect(face_handle_t const &f,
                                                      Vector const &v) const {
  Vector reflected{v};
  // Reverse the component corresponding to which plane the face lies in
  Plane p{face_to_plane(f.as_id())};
  switch (p) {
  case YZ: // X-component
    reflected[0] = -1.0 * reflected[0];
    break;
  case XZ: // Y-component
    reflected[1] = -1.0 * reflected[1];
    break;
  case XY: // Z-component
    reflected[2] = -1.0 * reflected[2];
    break;
  case NUM_PLANES:
    /* Wow should not be here*/
    // Insist(p < NUM_PLANES);
    printf("%s:%i reflect FAIL: unusable Plane!\n", __FUNCTION__, __LINE__);
    break;
  }
  return reflected;
}

namespace {
/**\brief The set of normals to faces, in order of Face_Name */
constexpr Vector face_normals[6] = {
    Vector{1.0, 0.0, 0.0},  Vector{0.0, 1.0, 0.0},  Vector{0.0, 0.0, 1.0},
    Vector{-1.0, 0.0, 0.0}, Vector{0.0, -1.0, 0.0}, Vector{0.0, 0.0, -1.0}};

} // namespace

inline Cartesian_Mesh::Face_Name
Cartesian_Mesh::face_to_face_name(cell_handle_t const &c,
                                  face_handle_t const &f) const {
  // Need to know direction, and whether it is high or low
  Plane p{face_to_plane(f.as_id())};
  // return the corresponding normal
  cell_handle_t possible_cell{face_to_cell(f)};
  bool is_high_face(possible_cell.as_id() != c.as_id());
  Face_Name fname{LOW_X};
  if (is_high_face) {
    switch (p) {
    case YZ:
      fname = HIGH_X;
      break;
    case XZ:
      fname = HIGH_Y;
      break;
    case XY:
      fname = HIGH_Z;
      break;
    case NUM_PLANES: // some hideous error here
      break;
    }
  } else {
    switch (p) {
    case YZ:
      fname = LOW_X;
      break;
    case XZ:
      fname = LOW_Y;
      break;
    case XY:
      fname = LOW_Z;
      break;
    case NUM_PLANES: // some hideous error here
      break;
    }
  }
  return fname;
}

inline std::vector<Cartesian_Mesh::face_handle_t>
Cartesian_Mesh::boundary_faces(Face_Name const f_n) const {
  std::vector<Cartesian_Mesh::face_handle_t> b_faces;
  index_t xmin = (f_n == HIGH_X) ? get_num_x() - 1u : 0u;
  index_t ymin = (f_n == HIGH_Y) ? get_num_y() - 1u : 0u;
  index_t zmin = (f_n == HIGH_Z) ? get_num_z() - 1u : 0u;
  index_t xmax = (f_n == LOW_X) ? 1u : get_num_x();
  index_t ymax = (f_n == LOW_Y) ? 1u : get_num_y();
  index_t zmax = (f_n == LOW_Z) ? 1u : get_num_z();
  // switch(case f_n) {
  //   case LOW_X:
  //     xmin = 0;
  //     xmax = 1;
  //     ymin = 0;
  //     ymax = get_num_y();
  //     zmin = 0;
  //     zmax = get_num_z();
  //     break;
  //   case LOW_Y:
  //     xmin = 0;
  //     xmax = get_num_x();
  //     ymin = 0;
  //     ymax = 1;
  //     zmin = 0;
  //     zmax = get_num_z();
  //     break;
  //   case LOW_Z:
  //     xmin = 0;
  //     xmax = get_num_x();
  //     ymin = 0;
  //     ymax = get_num_y();
  //     zmin = 0;
  //     zmax = 1;
  //     break;
  //   case HIGH_X:
  //     xmin = get_num_x() - 1;
  //     xmax = get_num_x();
  //     ymin = 0;
  //     ymax = get_num_y();
  //     zmin = 0;
  //     zmax = get_num_z();
  //     break;
  //   case HIGH_Y:
  //     xmin = 0;
  //     xmax = get_num_x();
  //     ymin = get_num_y() - 1;
  //     ymax = get_num_y();
  //     zmin = 0;
  //     zmax = get_num_z();
  //     break;
  //   case HIGH_Z:
  //     xmin = 0;
  //     xmax = get_num_x();
  //     ymin = 0;
  //     ymax = get_num_y();
  //     zmin = get_num_z() - 1;
  //     zmax = get_num_z();
  //     break;
  // }
  for (index_t iz = zmin; iz < zmax; ++iz) {
    for (index_t iy = ymin; iy < ymax; ++iy) {
      for (index_t ix = xmin; ix < xmax; ++ix) {
        cell_handle_t c{make_cell(ix, iy, iz)};
        face_handle_t f{cell_to_face(c, f_n)};
        b_faces.push_back(f);
      }
    }
  }
  return b_faces;
} // boundary_faces

inline Cartesian_Mesh::Vector
Cartesian_Mesh::get_normal(cell_handle_t const &c,
                           face_handle_t const &f) const {
  // get face name for this face in this cell
  Face_Name fname{face_to_face_name(c, f)};
  // Just compute the plane. Then use
  return face_normals[fname];
} // get_normal

inline Cartesian_Mesh::Vector Cartesian_Mesh::get_normal(Face_Name f) {
  return face_normals[f];
} // get_normal

template <typename RNG_T>
inline Cartesian_Mesh::Vector
Cartesian_Mesh::sample_position(RNG_T &rng, cell_handle_t const &c) const {
  auto [ix, iy, iz] = linear_to_cartesian(c.as_id(), nx_, ny_);
  geom_t delta_x{dx_ * rng.random()};
  geom_t delta_y{dy_ * rng.random()};
  geom_t delta_z{dz_ * rng.random()};
  return {xmin_ + ix * dx_ + delta_x, ymin_ + iy * dy_ + delta_y,
          zmin_ + iz * dz_ + delta_z};
} // sample_position

template <typename RNG_T>
inline Cartesian_Mesh::Vector
Cartesian_Mesh::sample_direction_isotropic(RNG_T &rng) {
  using nut::constants::pi;
  geom_t const ctheta = geom_t(2) * rng.random() - geom_t(1);
  geom_t const phi = geom_t(2) * pi * rng.random();
  geom_t const stheta = std::sqrt(1 - ctheta * ctheta);
  geom_t const cphi = std::cos(phi);
  geom_t const sphi = std::sin(phi);
  Vector v{{stheta * cphi, stheta * sphi, ctheta}};
  return v;
}

inline Cartesian_Mesh::cell_handle_t
Cartesian_Mesh::cell_across(cell_handle_t const &c, face_handle_t const &f,
                            Point const & /*p*/) const {
  auto [ix, iy, iz] = linear_to_cartesian(c.as_id(), nx_, ny_);
  switch (f.get_index()) {
  case LOW_X:
    ix = (ix > 0) ? ix - 1 : void_cell_idx;
    break;
  case HIGH_X:
    ix = (ix < nx_ - 1) ? ix + 1 : void_cell_idx;
    break;
  case LOW_Y:
    iy = (iy > 0) ? iy - 1 : void_cell_idx;
    break;
  case HIGH_Y:
    iy = (iy < ny_ - 1) ? iy + 1 : void_cell_idx;
    break;
  case LOW_Z:
    iz = (iz > 0) ? iz - 1 : void_cell_idx;
    break;
  case HIGH_Z:
    iz = (iz < nz_ - 1) ? iz + 1 : void_cell_idx;
    break;
    // default:
    // pithy error message here
  } // switch(face)
  index_t idx{};
  if (ix == void_cell_idx || iy == void_cell_idx || iz == void_cell_idx) {
    idx = void_cell_idx;
  } else {
    idx = cartesian_to_linear(ix, iy, iz, nx_, ny_);
  }
  return cell_handle_t{idx};
} // cell_across

inline Cartesian_Mesh::Geom_T
Cartesian_Mesh::volume(cell_handle_t const & /*c*/) const {
  return dx_ * dy_ * dz_;
}

/**\brief Find the cell in which the Point lies.
 *
 * Note that this will always return *something*. It will return junk if you
 * provide a point that is not in the mesh. You might want to check the
 * returned cell value with in_cell() (though you may get consistent junk).
 * So you may really want to check valid_cell.
 */
inline Cartesian_Mesh::cell_handle_t
Cartesian_Mesh::find_cell(Point const &p) const {
  index_t ix = std::floor((p[0] - xmin_) / dx_);
  index_t iy = std::floor((p[1] - ymin_) / dy_);
  index_t iz = std::floor((p[2] - zmin_) / dz_);
  cell_handle_t const c(cartesian_to_linear(ix, iy, iz, nx_, ny_));
  // Require(in_cell(p, c));
  return c;
} // find_cell

inline bool Cartesian_Mesh::valid_cell(cell_handle_t const &c) const {
  auto [ix, iy, iz] = linear_to_cartesian(c.as_id(), nx_, ny_);
  bool ok = ix < nx_ && iy < ny_ && iz < nz_;
  return ok;
}

/**\brief Compute the x-extents {low, high} for a cell. This is where the
 * actual encoding of the cell index is important--we treat the linear index
 * as encoding position information. */
inline Cartesian_Mesh::Extents_T
Cartesian_Mesh::get_x_extents(cell_handle_t const &c) const {
  auto [ix, iy, iz] = linear_to_cartesian(c.as_id(), nx_, ny_);
  geom_t x_lo = ix * dx_ + xmin_;
  geom_t x_hi = (ix + 1) * dx_ + xmin_;
  return {x_lo, x_hi};
}

/**\brief Compute the y-extents {low, high} for a cell. */
inline Cartesian_Mesh::Extents_T
Cartesian_Mesh::get_y_extents(cell_handle_t const &c) const {
  auto [ix, iy, iz] = linear_to_cartesian(c.as_id(), nx_, ny_);
  geom_t y_lo = iy * dy_ + ymin_;
  geom_t y_hi = (iy + 1) * dy_ + ymin_;
  return {y_lo, y_hi};
}

/**\brief Compute the z-extents {low, high} for a cell. */
inline Cartesian_Mesh::Extents_T
Cartesian_Mesh::get_z_extents(cell_handle_t const &c) const {
  auto [ix, iy, iz] = linear_to_cartesian(c.as_id(), nx_, ny_);
  geom_t z_lo = iz * dz_ + zmin_;
  geom_t z_hi = (iz + 1) * dz_ + zmin_;
  return {z_lo, z_hi};
}

inline geom_t Cartesian_Mesh::compute_distance(geom_t d_cos, geom_t coord,
                                               geom_t face_coord_lo,
                                               geom_t face_coord_hi) {
  using nut::constants::Max_Dbl;
  return d_cos != 0.0 ? (d_cos < 0.0 ? ((face_coord_lo - coord) / d_cos)
                                     : ((face_coord_hi - coord) / d_cos))
                      : Max_Dbl;
}

inline Cartesian_Mesh::intersect_t
Cartesian_Mesh::intersection(Ray const &r, cell_handle_t const &c) const {
  using nut::constants::Max_Dbl;
  Vector const &pos = r.position();
  // Require( in_cell(pos,c));
  Vector const &dir = r.direction();
  // Require(soft_equiv(dir.norm(),1.0));
  Face_Name f{LOW_X};
  geom_t dist = nut::constants::Max_Dbl;
  auto [x_lo, x_hi] = this->get_x_extents(c);
  auto [y_lo, y_hi] = this->get_y_extents(c);
  auto [z_lo, z_hi] = this->get_z_extents(c);
  const geom_t dx = compute_distance(dir[0], pos[0], x_lo, x_hi);
  const geom_t dy = compute_distance(dir[1], pos[1], y_lo, y_hi);
  const geom_t dz = compute_distance(dir[2], pos[2], z_lo, z_hi);
  if (dx <= dy) {
    dist = dx, f = dir[0] < 0.0 ? LOW_X : HIGH_X;
  } else {
    dist = dy, f = dir[1] < 0.0 ? LOW_Y : HIGH_Y;
  }
  if (dz < dist) {
    dist = dz, f = dir[2] < 0.0 ? LOW_Z : HIGH_Z;
  }
  // Check(dist < Max_Dbl);
  // Check(f < 6);

  /*Now , need to convert the face to an absolute face index*/
  return std::make_pair(Cartesian_Face(f), dist);
} // intersection

inline bool Cartesian_Mesh::in_cell(Vector const &p,
                                    cell_handle_t const &c) const {
  auto const [x_lo, x_hi] = get_x_extents(c);
  auto const [y_lo, y_hi] = get_y_extents(c);
  auto const [z_lo, z_hi] = get_z_extents(c);
  bool const x_ok = p[0] >= x_lo && p[0] < x_hi;
  bool const y_ok = p[1] >= y_lo && p[1] < y_hi;
  bool const z_ok = p[2] >= z_lo && p[2] < z_hi;
  return x_ok && y_ok && z_ok;
}

inline index_t Cartesian_Mesh::num_cells() const { return nx_ * ny_ * nz_; }

inline Cartesian_Mesh::cell_handle_t
Cartesian_Mesh::make_cell(index_t ix, index_t iy, index_t iz) const {
  return cell_handle_t{cartesian_to_linear(ix, iy, iz, nx_, ny_)};
}

// private methods

/**\brief Compute a Cartesian_Face index from a cell_handle_t index and one of
 * the six Face Names. Note: must be called with a mesh-valid cell; cannot be
 * called with +1 in any domension, for instance.
 */
inline Cartesian_Face Cartesian_Mesh::cell_to_face(cell_handle_t const &c,
                                                   Face_Name const fn) const {
  auto [ix, iy, iz] = linear_to_cartesian(c.as_id(), nx_, ny_);
  Plane p{to_plane(to_dir(fn))};
  if (is_high(fn)) {
    switch (p) {
    case YZ:
      ix++;
      break;
    case XZ:
      iy++;
      break;
    case XY:
      iz++;
      break;
    case NUM_PLANES:
      break;
      /* Shoule not be here!! */
    }
  }
  return cell_to_face(ix, iy, iz, p);
}

/**\brief The core cell_to_face computation.
 * Note: can be called with ix (iy, iz) greater than nx (ny, nz resp.) by one
 * and it will generally do the right thing.
 */
inline Cartesian_Face Cartesian_Mesh::cell_to_face(index_t ix, index_t iy,
                                                   index_t iz,
                                                   Plane const p) const {
  index_t base = 0;
  uint32_t ip = static_cast<uint32_t>(p);
  for (uint32_t i = 0; i < ip; ++i) {
    base += n_faces[i];
  }
  /* To get correct stride, if it's an X ("YZ") face, need to increment
   * nx; ditto for ny if it's a Y ("XZ") face */
  size_t const eff_nx = (p == YZ) ? nx_ + 1 : nx_;
  size_t const eff_ny = (p == XZ) ? ny_ + 1 : ny_;
  index_t arr_idx = cartesian_to_linear(ix, iy, iz, eff_nx, eff_ny);
  // Require(arr_idx < n_faces[p]);
  index_t the_index = base + arr_idx;
  // Require(the_index < base + n_faces[p]); // sort of restates the above,
  // eh?
  return Cartesian_Face(the_index);
} // cell_to_face

/**\brief Given a face index, to which plane does it belong? */
inline Plane Cartesian_Mesh::face_to_plane(index_t const face_idx) const {
  Plane plane{YZ}; // X face
  if (face_idx >= n_faces[1] + n_faces[0]) {
    plane = XY; // Z face
  } else if (face_idx >= n_faces[0]) {
    plane = XZ; // Y face
  }
  return plane;
}

/**\brief Given a face index, which */
inline index_t
Cartesian_Mesh::face_to_plane_offset(index_t const face_idx) const {
  // Find the greatest base address that can be subtracted off and this thing
  // still makes sense.
  index_t offset{0ull};
  if (face_idx >= n_faces[0]) {
    offset += n_faces[0];
  }
  if (face_idx >= n_faces[1] + n_faces[0]) {
    offset += n_faces[1];
  }
  return offset;
}

inline Cartesian_Mesh::cell_handle_t
Cartesian_Mesh::face_to_cell(Cartesian_Face const &f) const {
  index_t const full_idx = f.as_id();
  index_t const plane_offset{face_to_plane_offset(full_idx)};
  index_t const arr_idx{full_idx - plane_offset};
  cell_handle_t c{arr_idx};
  return c;
}

inline bool Cartesian_Mesh::is_boundary(face_handle_t f) const {
  /* Strategy: convert to cell, using the "face is always lower bound"
   * approach. Then, if the cell index in the corresponding dimension
   * is a boundary cell, it's a boundary face. */
  cell_handle_t const c{face_to_cell(f)};
  Plane p{face_to_plane(f.as_id())};
  index_t const eff_nx = (p == YZ) ? nx_ + 1 : nx_;
  index_t const eff_ny = (p == XZ) ? ny_ + 1 : ny_;
  index_t ix, iy, iz; // Intel doesn't do structured init til I18
  std::tie(ix, iy, iz) = linear_to_cartesian(c.as_id(), eff_nx, eff_ny);
  index_t relevant_idx = (p == YZ) ? ix : (p == XZ) ? iy : iz;
  index_t relevant_len = (p == YZ) ? nx_ : (p == XZ) ? ny_ : nz_;
  bool const is_bdy_idx = relevant_idx == 0 || relevant_idx == relevant_len;
  return is_bdy_idx;
} // is_boundary(face_handle_t )

inline Vertex Cartesian_Mesh::cell_to_vertex(index_t ix, index_t iy,
                                             index_t iz) const {
  return Vertex(cartesian_to_linear(ix, iy, iz, (nx_ + 1), (ny_ + 1)));
}

inline std::vector<Cartesian_Face>
Cartesian_Mesh::make_faces(index_t ix, index_t iy, index_t iz) const {
  /* Get the faces for the cell, convert them to absolute memory offset */
  std::vector<Cartesian_Face> fs{
      // order is low-x, low-y, low-z...
      cell_to_face(ix, iy, iz, YZ), cell_to_face(ix, iy, iz, XZ),
      cell_to_face(ix, iy, iz, XY),
      // ...then high-X, high-Y, high-Z
      cell_to_face(ix + 1, iy, iz, YZ), cell_to_face(ix, iy + 1, iz, XZ),
      cell_to_face(ix, iy, iz + 1, XY)};
  return fs;
} // make_faces

inline Cartesian_Mesh::vertex_list_t
Cartesian_Mesh::make_vertices(index_t ix, index_t iy, index_t iz) const {
  /* Get the faces for the cell, convert them to absolute memory offset */
  vertex_list_t vs = {cell_to_vertex(ix, iy, iz),
                      cell_to_vertex(ix + 1, iy, iz),
                      cell_to_vertex(ix, iy + 1, iz),
                      cell_to_vertex(ix + 1, iy + 1, iz),
                      cell_to_vertex(ix, iy, iz + 1),
                      cell_to_vertex(ix + 1, iy, iz + 1),
                      cell_to_vertex(ix, iy + 1, iz + 1),
                      cell_to_vertex(ix + 1, iy + 1, iz + 1)};
  return vs;
} // make_vertices

inline void Cartesian_Mesh::make_cartesian_mesh() {
  // set up the cells, faces, and vertices
  for (index_t iz = 0; iz < nz_; ++iz) {
    for (index_t iy = 0; iy < ny_; ++iy) {
      for (index_t ix = 0; ix < nx_; ++ix) {
        index_t cidx = cartesian_to_linear(ix, iy, iz, nx_, ny_);
        cell_handle_t c(cidx);
        std::vector<Cartesian_Face> fs(this->make_faces(ix, iy, iz));
        vertex_list_t vs(this->make_vertices(ix, iy, iz));
        this->add_cell(c, fs, vs);
      }
    }
  }
  return;
} // make_cartesian_mesh

// ------------------------------- iterator class ------------------------------
struct Cartesian_Mesh::boundary_face_iterator
    : public std::iterator<std::input_iterator_tag,
                           Cartesian_Mesh::face_handle_t> {
  using value_t = Cartesian_Mesh::face_handle_t;

  operator bool() const { return !at_end(); }

  bool at_end() const { return ix == xmax && iy == ymax && iz == zmax; }

  boundary_face_iterator &operator++() {
    ix++;
    if (ix == xmax) {
      iy++;
      if (iy == ymax) {
        iz++;
        if (at_end()) {
          return *this;
        }
        iy = ymin;
      }
      ix = xmin;
    }
    return *this;
  } // operator++

  boundary_face_iterator operator++(int) {
    boundary_face_iterator retval = *this;
    ++(*this);
    return retval;
  }

  bool operator==(boundary_face_iterator const &other) const {
    return ix == other.ix && iy == other.iy && iz == other.iz;
  }

  bool operator!=(boundary_face_iterator const &other) const {
    return !(*this == other);
  }

  value_t operator*() const {
    cell_handle_t c{m.make_cell(ix, iy, iz)};
    face_handle_t f{m.cell_to_face(c, f_n)};
    return f;
  }

  boundary_face_iterator(index_t ix_, index_t iy_, index_t iz_, index_t xmin_,
                         index_t xmax_, index_t ymin_, index_t ymax_,
                         index_t zmin_, index_t zmax_, Face_Name f_n_,
                         Cartesian_Mesh const &m_)
      : ix(ix_), iy(iy_), iz(iz_), xmin(xmin_), xmax(xmax_), ymin(ymin_),
        ymax(ymax_), zmin(zmin_), zmax(zmax_), f_n(f_n_), m(m_) {}

  // private:
  index_t ix;
  index_t iy;
  index_t iz;
  index_t xmin;
  index_t xmax;
  index_t ymin;
  index_t ymax;
  index_t zmin;
  index_t zmax;
  Face_Name f_n;
  Cartesian_Mesh const &m;
}; // boundary_face_iterator

inline Cartesian_Mesh::boundary_face_iterator
Cartesian_Mesh::boundary_faces_begin(Face_Name f) const {
  index_t xmin = (f == HIGH_X) ? get_num_x() - 1u : 0u;
  index_t ymin = (f == HIGH_Y) ? get_num_y() - 1u : 0u;
  index_t zmin = (f == HIGH_Z) ? get_num_z() - 1u : 0u;
  index_t xmax = (f == LOW_X) ? 1u : get_num_x();
  index_t ymax = (f == LOW_Y) ? 1u : get_num_y();
  index_t zmax = (f == LOW_Z) ? 1u : get_num_z();
  return boundary_face_iterator(xmin, ymin, zmin, xmin, xmax, ymin, ymax, zmin,
                                zmax, f, *this);
}

inline Cartesian_Mesh::boundary_face_iterator
Cartesian_Mesh::boundary_faces_end(Face_Name f) const {
  index_t xmin = (f == HIGH_X) ? get_num_x() - 1u : 0u;
  index_t ymin = (f == HIGH_Y) ? get_num_y() - 1u : 0u;
  index_t zmin = (f == HIGH_Z) ? get_num_z() - 1u : 0u;
  index_t xmax = (f == LOW_X) ? 1u : get_num_x();
  index_t ymax = (f == LOW_Y) ? 1u : get_num_y();
  index_t zmax = (f == LOW_Z) ? 1u : get_num_z();
  return boundary_face_iterator(xmax, ymax, zmax, xmin, xmax, ymin, ymax, zmin,
                                zmax, f, *this);
}

} // namespace nut_mesh

// End of file
