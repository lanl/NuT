// Mesh.hh
// T. M. Kelley
// Jan 11, 2011
// Header for Mesh
// (c) Copyright 2011 LANSLLC all rights reserved.

#ifndef MESH3DCar_H
#define MESH3DCar_H

#include "Assert.hh"
#include "constants.hh"
#include "lorentz.hh"
#include "soft_equiv.hh"
#include "types.hh"
#include "utilities.hh"
#include "utilities_io.hh"
#include <cmath>
#include <ostream>
#include <string>
#include <vector>

namespace nut {
/** \brief Index 3D things with a Cartesian tuple. */
template <typename cell_t>
struct ijk_t {
  cell_t i;
  cell_t j;
  cell_t k;
  ijk_t(cell_t i_, cell_t const j_, cell_t const k_) : i(i_), j(j_), k(k_) {}
};

/*!\brief Uniform 3D Cartesian mesh.
 * \tparam <cell_t> {cell index type}
 * \tparam <boundary_t> {geometry (numerical) type}
 * \tparam <bdy_descriptor_t> {boundary descriptor type}
 */
template <typename cell_t,
          typename geometry_t,
          typename bdy_descriptor_t = nut::bdy_types::descriptor>
struct Cartesian_3D {
public:
  // types and constants
  static const size_t dim = 3;

  using geom_t = geometry_t;
  using gcr = geom_t const &;
  using Vector = Vec_T<geom_t, dim>;

  using bdy_desc_t = bdy_descriptor_t;
  using vb = std::vector<geom_t>;
  using vec_bdy_descs = std::vector<bdy_desc_t>;
  using extents_t = std::pair<geom_t, geom_t>;
  using ijk_T = ijk_t<cell_t>;

  static const cell_t vac_cell = cell_t(-1);
  static constexpr geom_t huge = 1e300;

  static const cell_t null_cell_{0xFFFFFFFF};

  static cell_t null_cell() { return null_cell_; }

  /**\brief Position & direction */
  struct coord_t {
    vec_t<3> x;
    vec_t<3> omega;
  };  // coord_t

  enum face_t {
    low_x = 0, /* y-z plane */
    low_y,     /* x-z plane */
    low_z,     /* x-y plane */
    high_x,    /* y-z plane */
    high_y,    /* x-z plane */
    high_z     /* x-y plane */
  };

  using Face = face_t;

  enum dir_t { X = 0, Y, Z };

  /** \brief The answer to a distance to boundary calculation. */
  struct d_to_b_t {
    geom_t d;    /**< distance to boundary */
    face_t face; /**< which face will be intersected */
    d_to_b_t(geom_t const d_, face_t const f_) : d(d_), face(f_) {}
  };

  using Ray = coord_t;

  using Intersection = d_to_b_t;

  static geom_t get_distance(Intersection const & i) { return i.d; }

  static Face get_face(Intersection const & i) { return i.face; }

public:
  // Interface

  // ctor
  Cartesian_3D(cell_t const nx,
               cell_t const ny,
               cell_t const nz,
               geom_t const dx,
               geom_t const dy,
               geom_t const dz,
               vec_bdy_descs const bds);

  cell_t n_cells() const { return m_ncells; }

  /*!\brief volume of cell. */
  geom_t volume(cell_t const /* ignored! */) const
  {
    geom_t vol = m_dx * m_dy * m_dz;
    return vol;
  }  // volume

  /** Is a cell a boundary for the relevant face? */
  inline bool isBoundary(cell_t const cidx, face_t const face) const;

  /*!\brief which cell is across a face. Computed.
   *
   * \param cell: 0 < cell <= n_cells
   * \param face: enum instances
   */
  inline cell_t cell_across(cell_t const cell, face_t face) const;

  /** \brief Sample a position in a cell. */
  template <typename RNG_T>
  inline Vector sample_position(RNG_T & rng, cell_t const cell) const;

  /** Sample a direction uniformly on the unit sphere. */
  template <typename RNG_T>
  static Vector sample_direction(RNG_T & rng);

  /*!\brief get the lower and upper bounds of a cell in the direction indicated
   * by the face. For example, passing low or high x gives the x direction. *
   * 0 < cell <= n_cells                              */
  extents_t cell_extents(cell_t const cell, dir_t const d) const;

  /*!\brief calculate new coordinates and new direction cosines at a given
   *        distance along old direction cosines.  */
  coord_t new_coordinate(vec_t<3> const & oldx,
                         vec_t<3> const & oldomega,
                         geom_t const distance) const;

  inline Intersection intersection(Ray const & r, cell_t const c) const;

  /** Distance-to-boundary. This is possibly not the most elegant
   * or efficient or vectorizable implementation. */
  inline d_to_b_t distance_to_bdy(coord_t const & crd, cell_t const cell) const;

  inline d_to_b_t distance_to_bdy(vec_t<3> const & x,
                                  vec_t<3> const & omega,
                                  cell_t const cell) const;

  /** This stores boundary conditions for a rectangular solid mesh.
   * BCs are loaded into the vector of
   * boundary descriptors in the following order:
   * low-x plane (y-z)
   * high-x plane
   * low-y plane (x-z)
   * high-y plane
   * low-z plane (x-y)
   * high-z plane
   *
   * In each plane, b.c.'s are loaded with the lower coordinate
   * running first (contiguous). To recover the boundary condition
   * for a cell, first convert to i,j,k, then use the correct offset
   * and the correct two indices.
   */
  inline bdy_desc_t getBC(face_t const face, cell_t const cidx) const;

  /*!\brief type of cell boundary */
  inline bdy_desc_t get_bdy_type(cell_t const c, cell_t const face) const;

  /** Transform energy and momentum to the frame that is moving with
   * velocity v in the lab frame. v is the velocity of the co-moving
   * frame in the lab frame. */
  static EandOmega<3> LT_to_comoving(vec_t<3> const & v_lab,
                                     geom_t const & e_lab,
                                     vec_t<3> const & omega_lab)
  {
    return spec_3D_Cartesian::LT_to_comoving(v_lab, e_lab, omega_lab);
  }  // LT_to_comoving

  /** Transform energy and momentum to the lab frame from the
      co-moving frame. Note that v is the velocity of the co-moving
      frame in the lab frame (not the velocity of the lab frame in c-m.
      This accomodates storing all the material velocities in the
      lab frame. */
  static EandOmega<3> LT_to_lab(vec_t<3> const & v_lab,
                                geom_t const & e_com,
                                vec_t<3> const & omega_com)
  {
    return spec_3D_Cartesian::LT_to_comoving(v_lab, e_com, omega_com);
  }  // LT_to_comoving

private:
  // state
  cell_t const m_nx;
  cell_t const m_ny;
  cell_t const m_nz;
  geom_t const m_dx;
  geom_t const m_dy;
  geom_t const m_dz;

  cell_t const m_ncells;
  vec_bdy_descs m_bdy_descs;

  // implementation
private:
  /** Get the grid spacing in the specified direction*/
  inline geom_t getDeltaCoord(dir_t const d) const
  {
    geom_t res(0);
    switch(d) {
      case X: res = m_dx; break;
      case Y: res = m_dy; break;
      case Z: res = m_dz; break;
    }
    return res;
  }

  static inline ijk_T updateIJK(ijk_T const & old,
                                cell_t const newIdx,
                                face_t const f)
  {
    ijk_T out(old.i, old.j, old.k);
    switch(f % 3) {
      case 0: out.i = newIdx; break;
      case 1: out.j = newIdx; break;
      case 2: out.k = newIdx; break;
      default: Require(false, "Mesh3DCar::updateIJK: face out of range!");
    }
    return out;
  }

  /** \brief True if face is on the high side. */
  static inline bool isHigh(face_t const face) { return face / 3; }

  static inline cell_t getIJKComp(ijk_T const & ijk, face_t const face)
  {
    return getIJKComp(ijk, dir_t(face % 3));
  }
  /** \brief Get copy of Cartesian component. */
  static inline cell_t getIJKComp(ijk_T const & ijk, dir_t const d)
  {
    cell_t result(0);
    switch(d) {
      case X: result = ijk.i; break;
      case Y: result = ijk.j; break;
      case Z: result = ijk.k; break;
      default: Require(false, "getIJKComp: direction out of range???");
    }
    return result;
  }

  inline cell_t getDimSize(face_t const face) const
  {
    if(low_x == face || high_x == face) return m_nx;
    if(low_y == face || high_y == face) return m_ny;
    return m_nz;
  }

  /** \brief Convert a lexical cell index to (i,j,k) */
  inline ijk_T mkIJK(cell_t const idx) const
  {
    ijk_T ijk((idx % m_nx), (idx / m_nx % m_ny), (idx / m_nx / m_ny));
    return ijk;
  }
  /** \brief Convert Cartesian index to lexical. */
  inline cell_t toIndex(ijk_T const ijk) const
  {
    return ijk.i + ijk.j * m_nx + ijk.k * m_nx * m_ny;
  }

  static void printIJK(ijk_T const & ijk)
  {
    std::cout << ijk.i << "," << ijk.j << "," << ijk.k << "\n";
  }

  static inline cell_t getBCIndex(cell_t const i1,
                                  cell_t const i2,
                                  cell_t const n1)
  {
    return i1 + n1 * i2;
  }

public:
  /* This is forced to be public, b/c external functions use it, and it has
   * to be a member because of the types. */
  static inline cell_t getBCOffset(face_t const face,
                                   cell_t const nx,
                                   cell_t const ny,
                                   cell_t const nz)
  {
    cell_t offset(0);
    /* The offset for a given face is the sum of all the previous
     * offsets, so yes, fall through cases. That is, read this
     * from bottom to top, stopping at the face of interest. */
    switch(face) {
      case high_z: offset += nx * ny;  // + nx * ny for low-z
      case low_z: offset += nx * nz;   // + nx * nz for high-y
      case high_y: offset += nx * nz;  // + nx * nz for low-y
      case low_y: offset += ny * nz;   // + ny * nz offsets for high-X
      case high_x: offset += ny * nz;  // there were ny * nz offsets for low-X
      case low_x: offset += 0;         // offsets start at 0 with lox-x
    }
    return offset;
  }  // getBCOffset
};   // struct Cartesian_3D

// template <typename cell_t, typename geometry_t, typename bdy_descriptor_t>
// const cell_t Cartesian_3D<cell_t, geometry_t, bdy_descriptor_t>::null_cell_ =
//     0xFFFFFFFF;

/** \brief Make a plane's worth of reflecting b.c.'s; helper for mkReflectBCs.
 */
template <typename cell_t,
          typename bdy_descriptor_t = nut::bdy_types::descriptor>
cell_t
mkPlaneBC(std::vector<bdy_descriptor_t> & bds,
          cell_t const n1,
          cell_t const n2,
          cell_t const offset,
          bdy_descriptor_t const bdy)
{
  for(cell_t i2 = 0; i2 < n2; ++i2) {
    for(cell_t i1 = 0; i1 < n1; ++i1) {
      cell_t const idx = i1 + n1 * i2 + offset;
      bds[idx] = bdy;
    }
  }
  return n1 * n2;
}

/** \brief Make a set of boundary conditions that are reflecting on every face.
 */
template <typename mesh_t,
          typename cell_t,
          typename bdy_descriptor_t = nut::bdy_types::descriptor>
void
mkReflectBCs(std::vector<bdy_descriptor_t> & bds,
             cell_t const nx,
             cell_t const ny,
             cell_t const nz)
{
  cell_t const n_bs = 2 * nx * ny + 2 * nx * nz + 2 * ny * nz;
  bdy_descriptor_t const reflect{nut::bdy_types::REFLECTIVE};
  if(n_bs != bds.size()) { bds.resize(n_bs, reflect); }
  cell_t offset(0);
  // low x bounds
  Require(offset == mesh_t::getBCOffset(mesh_t::low_x, nx, ny, nz),
          "mkReflectBCs: offset not correct 514");
  offset += mkPlaneBC(bds, ny, nz, offset, reflect);
  // high x bounds
  Require(offset == mesh_t::getBCOffset(mesh_t::high_x, nx, ny, nz),
          "mkReflectBCs: offset not correct 518");
  offset += mkPlaneBC(bds, ny, nz, offset, reflect);
  // low y bounds
  Require(offset == mesh_t::getBCOffset(mesh_t::low_y, nx, ny, nz),
          "mkReflectBCs: offset not correct 522");
  offset += mkPlaneBC(bds, nx, nz, offset, reflect);
  // high y bounds
  Require(offset == mesh_t::getBCOffset(mesh_t::high_y, nx, ny, nz),
          "mkReflectBCs: offset not correct 526");
  offset += mkPlaneBC(bds, nx, nz, offset, reflect);
  // low z bounds
  Require(offset == mesh_t::getBCOffset(mesh_t::low_z, nx, ny, nz),
          "mkReflectBCs: offset not correct 530");
  offset += mkPlaneBC(bds, nx, ny, offset, reflect);
  // high z bounds
  Require(offset == mesh_t::getBCOffset(mesh_t::high_z, nx, ny, nz),
          "mkReflectBCs: offset not correct 534");
  offset += mkPlaneBC(bds, nx, ny, offset, reflect);

  Require(offset == bds.size(),
          "mkReflectBCs: offset did not account for all entries");
  return;
}

template <typename mesh_t,
          typename cell_t,
          typename bdy_descriptor_t = nut::bdy_types::descriptor>
void
mkLowRefHighVacBounds(std::vector<nut::bdy_types::descriptor> & bds,
                      cell_t const nx,
                      cell_t const ny,
                      cell_t const nz)
{
  cell_t const n_bs = 2 * nx * ny + 2 * nx * nz + 2 * ny * nz;
  nut::bdy_types::descriptor ref(nut::bdy_types::REFLECTIVE);
  nut::bdy_types::descriptor vac(nut::bdy_types::VACUUM);
  nut::bdy_types::descriptor trn(nut::bdy_types::CELL);

  if(n_bs != bds.size()) { bds.resize(n_bs); }
  bds.assign(bds.size(), trn);
  cell_t offset(0);

  // low z bounds
  offset += nut::mkPlaneBC(bds, nx, ny, offset, ref);
  // high z bounds
  Require(offset == mesh_t::getBCOffset(mesh_t::high_z, nx, ny, nz),
          "offset not correct");
  offset += nut::mkPlaneBC(bds, nx, ny, offset, vac);
  // low y bounds
  Require(offset == mesh_t::getBCOffset(mesh_t::low_y, nx, ny, nz),
          "mkReflectBCs: offset not correct");
  offset += nut::mkPlaneBC(bds, nx, nz, offset, ref);
  // high y bounds
  Require(offset == mesh_t::getBCOffset(mesh_t::high_y, nx, ny, nz),
          "mkReflectBCs: offset not correct");
  offset += nut::mkPlaneBC(bds, nx, nz, offset, vac);
  // low x bounds
  Require(offset == mesh_t::getBCOffset(mesh_t::low_x, nx, ny, nz),
          "mkReflectBCs: offset not correct");
  offset += nut::mkPlaneBC(bds, ny, nz, offset, ref);
  // high x bounds
  Require(offset == mesh_t::getBCOffset(mesh_t::high_x, nx, ny, nz),
          "mkReflectBCs: offset not correct");
  offset += nut::mkPlaneBC(bds, ny, nz, offset, vac);
  return;
}  // mkLowRefHighVacBounds

// useful typedef
using Cartesian_3D_Mesh = Cartesian_3D<cell_t, geom_t>;

Cartesian_3D_Mesh
make_vacuum_bdy_c3()
{
  cell_t const nx = 3, ny = 5, nz = 7;
  geom_t const dx = 1.0, dy = 2.0, dz = 3.0;
  Cartesian_3D_Mesh::vec_bdy_descs bds{};
  Cartesian_3D_Mesh mesh(nx, ny, nz, dx, dy, dz, bds);
  return mesh;
}

}  // namespace nut

// introduce inline template definitions
#define ALLOWED_TO_INCLUDE_MESH3D_CAR_I_HH
#include "Mesh3DCar.i.hh"
#undef ALLOWED_TO_INCLUDE_MESH3D_CAR_I_HH

#endif

// version
// $Id$

// End of file
