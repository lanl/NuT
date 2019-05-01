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
  static const size_t dim = 3;

  typedef geometry_t geom_t;
  typedef geom_t const & gcr;
  // typedef Vec_T<geom_t,3> vec_t;

  typedef bdy_descriptor_t bdy_desc_t;
  typedef std::vector<geom_t> vb;
  typedef std::vector<bdy_desc_t> vbd;
  typedef std::pair<geom_t, geom_t> extents_t;
  typedef ijk_t<cell_t> ijk_T;

  static const cell_t vac_cell = cell_t(-1);
  static constexpr geom_t huge = 1e300;

  /**\brief Position & direction */
  struct coord_t {
    vec_t<3> x;
    vec_t<3> omega;
  };  // coord_t

  enum face_t {
    low_x = 0, /* y-z plane */
    high_x,    /* y-z plane */
    low_y,     /* x-z plane */
    high_y,    /* x-z plane */
    low_z,     /* x-y plane */
    high_z     /* x-y plane */
  };

  enum dir_t { X = 0, Y, Z };

  /** \brief The answer to a distance to boundary calculation. */
  struct d_to_b_t {
    geom_t d;    /**< distance to boundary */
    face_t face; /**< which face will be intersected */
    d_to_b_t(geom_t const d_, face_t const f_) : d(d_), face(f_) {}
  };

  /** \brief Convert a lexical cell index to (i,j,k) */
  inline ijk_T mkIJK(cell_t const idx) const
  {
    ijk_T ijk((idx % m_nx), (idx / m_nx % m_ny), (idx / m_nx / m_ny));
    return std::move(ijk);
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

  static inline cell_t getBCOffset(face_t const face,
                                   cell_t const nx,
                                   cell_t const ny,
                                   cell_t const nz)
  {
    cell_t offset(0);
    /* The offset for a given face is the sum of all the
     * previous offsets, so yes, fall through cases. */
    switch(face) {
      case high_x: offset += ny * nz;
      case low_x: offset += nx * nz;
      case high_y: offset += nx * nz;
      case low_y: offset += nx * ny;
      case high_z: offset += nx * ny;
      case low_z: offset += 0;
    }
    return offset;
  }

  /** This uses a cheap and dirty way of storing boundary conditions *
   * for a rectangular solid. BCs are loaded into the vector of *
   * boundary descriptors in the following order:
   * low-z plane (x-y)
   * high-z plane
   * low-y plane (x-z)
   * high-y plane
   * low-x plane (y-z)
   * high-x plane
   *
   * In each plane, b.c.'s are loaded with the lower coordinate
   * running first (contiguous). To recover the boundary condition
   * for a cell, first convert to i,j,k, then use the correct offset
   * and the correct two indices.
   */
  bdy_desc_t getBC(face_t const face, cell_t const cidx) const
  {
    cell_t const offset(getBCOffset(face, m_nx, m_ny, m_nz));
    ijk_T const ijk(mkIJK(cidx));
    cell_t const i1 = (face >= low_y) ? ijk.i : ijk.k;
    cell_t const i2 = (face <= high_y) ? ijk.k : ijk.k;
    cell_t const n1 = (face <= high_y) ? m_nx : m_ny;
    cell_t const index = offset + getBCIndex(i1, i2, n1);
    return m_bds[index];
  }

  // ctor
  Cartesian_3D(geom_t const dx,
               cell_t const nx,
               geom_t const dy,
               cell_t const ny,
               geom_t const dz,
               cell_t const nz,
               vbd const bds)
      : m_dx(dx),
        m_nx(nx),
        m_dy(dy),
        m_ny(ny),
        m_dz(dz),
        m_nz(nz),
        m_ncells(nx * ny * nz),
        m_bds(bds)
  {
    nut::Require(m_bds.size() == 2 * nx * ny + 2 * nx * nz + 2 * ny * nz,
                 "Cartesian_3D: Must have correct number of"
                 " boundary descriptors");
    auto gtg = [](geom_t const x, geom_t const e, const char * s) {
      return nut::GreaterThan(x, e, s);
    };
    auto gtc = [](cell_t const x, cell_t const e, const char * s) {
      return nut::GreaterThan(x, e, s);
    };
    gtg(dx, 0.0, "Cartesian_3D: x cell spacing must be > 0");
    gtg(dy, 0.0, "Cartesian_3D: y cell spacing must be > 0");
    gtg(dz, 0.0, "Cartesian_3D: z cell spacing must be > 0");
    gtc(nx, 0, "Cartesian_3D: nx must be > 0");
    gtc(ny, 0, "Cartesian_3D: ny must be > 0");
    gtc(nz, 0, "Cartesian_3D: nz must be > 0");
  }

  cell_t n_cells() const { return m_ncells; }

  /*!\brief volume of cell. */
  geom_t volume(cell_t const /* ignored! */) const
  {
    geom_t vol = m_dx * m_dy * m_dz;
    return vol;
  }  // volume
  /** \brief True if face is on the high side. */
  static inline bool isHigh(face_t const face) { return face % 2; }

  static inline cell_t getIJKComp(ijk_T const & ijk, face_t const face)
  {
    return getIJKComp(ijk, dir_t(face / 2));
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

  /** Is a cell a boundary for the relevant face? */
  inline bool isBoundary(cell_t const cidx, face_t const face) const
  {
    ijk_T const ijk(mkIJK(cidx));
    return isBoundary(ijk, face);
  }

  inline bool isBoundary(ijk_T const & ijk, face_t const face) const
  {
    bool isBound(false);
    cell_t const idx = getIJKComp(ijk, face);
    bool const ishigh(isHigh(face));
    if(0 == idx && !ishigh)
      isBound = true;
    else {
      cell_t ndir = getDimSize(face);
      if(ndir - 1 == idx && ishigh) isBound = true;
    }
    return isBound;
  }

  static inline ijk_T updateIJK(ijk_T const & old,
                                cell_t const newIdx,
                                face_t const f)
  {
    ijk_T out(old.i, old.j, old.k);
    switch(f / 2) {
      case 0: out.i = newIdx; break;
      case 1: out.j = newIdx; break;
      case 2: out.k = newIdx; break;
      default: Require(false, "Mesh3DCar::updateIJK: face out of range!");
    }
    return std::move(out);
  }

  /*!\brief which cell is across a face. Computed.
   *
   * \param cell: 0 < cell <= n_cells
   * \param face: enum instances
   */
  cell_t cell_across_face(cell_t const cell, face_t face) const
  {
    // convert to i,j,k
    ijk_T const ijkIn(mkIJK(cell));
    cell_t result(cell_t(0));
    // add/subtract 1 from appropriate index
    if(!isBoundary(cell, face)) {
      cell_t inew = getIJKComp(ijkIn, face) + (isHigh(face) ? 1 : (-1));
      ijk_T updd = updateIJK(ijkIn, inew, face);
      result = toIndex(updd);
    }
    else {
      // otherwise, consult boundary condition
      bdy_desc_t const bd = this->getBC(face, cell);
      switch(bd) {
        case nut::bdy_types::R: result = cell; break;
        case nut::bdy_types::V: result = vac_cell; break;
        default:
          std::stringstream errstr;
          errstr << "Mesh3DCar::cell_across_face: "
                 << "unexpected boundary descriptor " << bd;
      }
    }
    return result;
  }  // cell_across_face

  /** \brief Sample a position in a cell. */
  template <typename RNG_T>
  inline coord_t sample_position(RNG_T & rng, cell_t const cell) const
  {
    coord_t coord;
    extents_t xs = cell_extents(cell, X);
    coord.x.v[0] = xs.first + rng.random() * (xs.second - xs.first);
    extents_t ys = cell_extents(cell, Y);
    coord.x.v[1] = ys.first + rng.random() * (ys.second - ys.first);
    extents_t zs = cell_extents(cell, Z);
    coord.x.v[2] = zs.first + rng.random() * (zs.second - zs.first);
    coord.omega = sample_direction(rng);
    return std::move(coord);
  }  // sample_position

  /** Sample a direction uniformly on the unit sphere. */
  template <typename RNG_T>
  static inline vec_t<3> sample_direction(RNG_T & rng)
  {
    geom_t const ctheta = geom_t(2) * rng.random() - geom_t(1);
    geom_t const phi = geom_t(2) * pi * rng.random();
    geom_t const stheta = std::sqrt(1 - ctheta * ctheta);
    geom_t const cphi = std::cos(phi);
    geom_t const sphi = std::sin(phi);
    vec_t<3> v;
    v.v = {{stheta * cphi, stheta * sphi, ctheta}};
    return v;
  }

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

  /*!\brief get the lower and upper bounds of a cell in the direction indicated
   * by the face. For example, passing low or high x gives the x direction. *
   * 0 < cell <= n_cells                              */
  extents_t cell_extents(cell_t const cell, dir_t const d) const
  {
    ijk_T ijk(mkIJK(cell));
    cell_t const idx = getIJKComp(ijk, d);
    cell_t const dc = getDeltaCoord(d);
    geom_t clow = idx * dc;
    geom_t chigh = (idx + 1) * dc;
    return extents_t(clow, chigh);
    // return extents_t(m_bdys[cell-1],m_bdys[cell]);
  }

  /*!\brief calculate new coordinates and new direction cosines at a given
   *        distance along old direction cosines.  */
  coord_t new_coordinate(coord_t const oldc, geom_t const distance) const
  {
    return this->new_coordinate(oldc.x, oldc.omega, distance);
  }

  coord_t new_coordinate(vec_t<3> const & oldx,
                         vec_t<3> const & oldomega,
                         geom_t const distance) const
  {
    coord_t newc;
    newc.x.v[0] = oldx.v[0] + distance * oldomega.v[0];
    newc.x.v[1] = oldx.v[1] + distance * oldomega.v[1];
    newc.x.v[2] = oldx.v[2] + distance * oldomega.v[2];
    newc.omega.v[0] = oldomega.v[0];
    newc.omega.v[1] = oldomega.v[1];
    newc.omega.v[2] = oldomega.v[2];
    return std::move(newc);
  }  // new_coordinate

  /** Distance-to-boundary. This is possibly not the most elegant
   * or efficient or vectorizable implementation. */
  d_to_b_t distance_to_bdy(coord_t const crd, cell_t const cell) const
  {
    return this->distance_to_bdy(crd.x, crd.omega, cell);
  }

  d_to_b_t distance_to_bdy(vec_t<3> const & x,
                           vec_t<3> const & omega,
                           cell_t const cell) const
  {
    extents_t const xs = cell_extents(cell, X);
    d_to_b_t const dx =
        omega.v[0] > 0
            ? d_to_b_t((xs.second - x.v[0]) / omega.v[0], high_x)
            : omega.v[0] < 0 ? d_to_b_t((xs.first - x.v[0]) / omega.v[0], low_x)
                             : d_to_b_t(huge, high_x);
    extents_t const ys = cell_extents(cell, Y);
    d_to_b_t const dy =
        omega.v[1] > 0
            ? d_to_b_t((ys.second - x.v[1]) / omega.v[1], high_y)
            : omega.v[1] < 0 ? d_to_b_t((ys.first - x.v[1]) / omega.v[1], low_y)
                             : d_to_b_t(huge, high_y);

    extents_t const zs = cell_extents(cell, Z);
    d_to_b_t const dz =
        omega.v[2] > 0
            ? d_to_b_t((zs.second - x.v[2]) / omega.v[2], high_z)
            : omega.v[2] < 0 ? d_to_b_t((zs.first - x.v[2]) / omega.v[2], low_z)
                             : d_to_b_t(huge, high_z);

    auto comp = [](d_to_b_t const * d1, d_to_b_t const * d2) {
      return (d1->d) < (d2->d) ? true : false;
    };

    d_to_b_t const * result = std::min(&dx, std::min(&dy, &dz, comp), comp);
    return std::move(*result);
  }

  geom_t const m_dx;
  cell_t const m_nx;
  geom_t const m_dy;
  cell_t const m_ny;
  geom_t const m_dz;
  cell_t const m_nz;

  cell_t const m_ncells;
  vbd m_bds;

  struct EAndOmega {
    vec_t<3> omega;
    geom_t e;
  };

  static geom_t inline gamma(vec_t<3> const v)
  {
    return 1.0 / std::sqrt(1 - dot(v, v) / (c * c));
  }

  /** Transform energy and momentum to the frame that is moving with
   * velocity v in the lab frame. v is the velocity of the co-moving
   * frame in the lab frame. */
  static EandOmega<3> inline LT_to_comoving(vec_t<3> const & v_lab,
                                            geom_t const & e_lab,
                                            vec_t<3> const & omega_lab)
  {
    EAndOmega res;
    geom_t const vdo = dot(v_lab, omega_lab);
    geom_t const gam = gamma(v_lab);
    geom_t const goc = gam / c;
    geom_t const fac = 1 - goc * vdo / (gam + 1);
    res.e = gam * e_lab * (1 - vdo / c);
    geom_t const eoec = e_lab / res.e;  // E/E_com
    res.omega.v[0] = eoec * (omega_lab.v[0] - goc * v_lab.v[0] * fac);
    res.omega.v[1] = eoec * (omega_lab.v[1] - goc * v_lab.v[1] * fac);
    res.omega.v[2] = eoec * (omega_lab.v[2] - goc * v_lab.v[2] * fac);
    return std::move(res);
  }  // LT_to_comoving

  /** Transform energy and momentum to the lab frame from the
      co-moving frame. Note that v is the velocity of the co-moving
      frame in the lab frame (not the velocity of the lab frame in c-m.
      This accomodates storing all the material velocities in the
      lab frame. */
  static EandOmega<3> inline LT_to_lab(vec_t<3> const & v_lab,
                                       geom_t const & e_com,
                                       vec_t<3> const & omega_com)
  {
    EandOmega<3> res;
    geom_t const vdo = dot(v_lab, omega_com);
    geom_t const gam = gamma(v_lab);
    geom_t const goc = gam / c;
    geom_t const fac = 1 + goc * vdo / (gam + 1);
    res.first = gam * e_com * (1 + vdo / c);
    geom_t const eoec = e_com / res.first;  // E/E_com
    res.second.v[0] = eoec * (omega_com.v[0] + goc * v_lab.v[0] * fac);
    res.second.v[1] = eoec * (omega_com.v[1] + goc * v_lab.v[1] * fac);
    res.second.v[2] = eoec * (omega_com.v[2] + goc * v_lab.v[2] * fac);
    return res;
  }  // LT_to_comoving

};  // Cartesian_3D

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
  bdy_descriptor_t const ref(nut::bdy_types::R);
  cell_t const n_bs = 2 * nx * ny + 2 * nx * nz + 2 * ny * nz;
  if(n_bs != bds.size()) { bds.resize(n_bs); }
  cell_t offset(0);
  // low z bounds
  offset += mkPlaneBC(bds, nx, ny, offset, ref);
  // high z bounds
  Require(offset == mesh_t::getBCOffset(mesh_t::high_z, nx, ny, nz),
          "offset not correct");
  offset += mkPlaneBC(bds, nx, ny, offset, ref);
  // low y bounds
  Require(offset == mesh_t::getBCOffset(mesh_t::low_y, nx, ny, nz),
          "mkReflectBCs: offset not correct");
  offset += mkPlaneBC(bds, nx, nz, offset, ref);
  // high y bounds
  Require(offset == mesh_t::getBCOffset(mesh_t::high_y, nx, ny, nz),
          "mkReflectBCs: offset not correct");
  offset += mkPlaneBC(bds, nx, nz, offset, ref);
  // low x bounds
  Require(offset == mesh_t::getBCOffset(mesh_t::low_x, nx, ny, nz),
          "mkReflectBCs: offset not correct");
  offset += mkPlaneBC(bds, ny, nz, offset, ref);
  // high x bounds
  Require(offset == mesh_t::getBCOffset(mesh_t::high_x, nx, ny, nz),
          "mkReflectBCs: offset not correct");
  offset += mkPlaneBC(bds, ny, nz, offset, ref);
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
  nut::bdy_types::descriptor ref(nut::bdy_types::R);
  nut::bdy_types::descriptor vac(nut::bdy_types::V);
  nut::bdy_types::descriptor trn(nut::bdy_types::T);

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

}  // namespace nut

#endif

// version
// $Id$

// End of file
