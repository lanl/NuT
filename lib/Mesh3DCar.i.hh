// Mesh3DCar.i.hh
// T. M. Kelley
// May 14, 2019
// (c) Copyright 2019 Triad National Security, all rights reserved

/* Implementations for Mesh3DCar methods */

#ifndef ALLOWED_TO_INCLUDE_MESH3D_CAR_I_HH
#error "Do not include Mesh3DCar.i.hh direectly: instead include Mesh3DCar.hh"
#endif

#pragma once

namespace nut {

template <typename cell_t, typename geometry_t, typename bdy_descriptor_t>
typename Cartesian_3D<cell_t, geometry_t, bdy_descriptor_t>::Vector
Cartesian_3D<cell_t, geometry_t, bdy_descriptor_t>::get_normal(Face const & f)
{
  Vector n;
  switch(f.id()) {
    case low_x: n = Vector{1.0, 0.0, 0.0}; break;
    case low_y: n = Vector{0.0, 1.0, 0.0}; break;
    case low_z: n = Vector{0.0, 0.0, 1.0}; break;
    case high_x: n = Vector{-1.0, 0.0, 0.0}; break;
    case high_y: n = Vector{0.0, -1.0, 0.0}; break;
    case high_z: n = Vector{0.0, 0.0, -1.0}; break;
  }
  return n;
}  // get_normal

template <typename cell_t, typename geometry_t, typename bdy_descriptor_t>
typename Cartesian_3D<cell_t, geometry_t, bdy_descriptor_t>::Vector
Cartesian_3D<cell_t, geometry_t, bdy_descriptor_t>::reflect(Face const & f,
                                                            Vector const & v)
{
  Vector normal{get_normal(f)};
  Vector reflected = v.reflect(normal);
  return reflected;
}

template <typename cell_t, typename geometry_t, typename bdy_descriptor_t>
bdy_descriptor_t
Cartesian_3D<cell_t, geometry_t, bdy_descriptor_t>::get_boundary_type(
    cell_t const c,
    cell_t const face) const
{
  cell_t const idx = make_idx(c, m_ncells);
  return m_bdy_descs[idx + face];
}

template <typename cell_t, typename geometry_t, typename bdy_descriptor_t>
typename Cartesian_3D<cell_t, geometry_t, bdy_descriptor_t>::Intersection
Cartesian_3D<cell_t, geometry_t, bdy_descriptor_t>::intersection(
    Ray const & r,
    cell_t const c) const
{
  vec_t<dim> & x{r.position()};
  vec_t<dim> & o{r.direction()};
  return this->distance_to_bdy(x, o, c);
}

template <typename cell_t, typename geometry_t, typename bdy_descriptor_t>
typename Cartesian_3D<cell_t, geometry_t, bdy_descriptor_t>::d_to_b_t
    Cartesian_3D<cell_t, geometry_t, bdy_descriptor_t>::distance_to_bdy(
        vec_t<3> const & x,
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
  return *result;
}

template <typename cell_t, typename geometry_t, typename bdy_descriptor_t>
typename Cartesian_3D<cell_t, geometry_t, bdy_descriptor_t>::d_to_b_t
Cartesian_3D<cell_t, geometry_t, bdy_descriptor_t>::distance_to_bdy(
    coord_t const & crd,
    cell_t const cell) const
{
  return this->distance_to_bdy(crd.x, crd.omega, cell);
}

template <typename cell_t, typename geometry_t, typename bdy_descriptor_t>
typename Cartesian_3D<cell_t, geometry_t, bdy_descriptor_t>::coord_t
    Cartesian_3D<cell_t, geometry_t, bdy_descriptor_t>::new_coordinate(
        vec_t<3> const & oldx,
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
  return newc;
}  // new_coordinate

template <typename cell_t, typename geometry_t, typename bdy_descriptor_t>
typename Cartesian_3D<cell_t, geometry_t, bdy_descriptor_t>::extents_t
Cartesian_3D<cell_t, geometry_t, bdy_descriptor_t>::cell_extents(
    cell_t const cell,
    dir_t const d) const
{
  ijk_T ijk(mkIJK(cell));
  cell_t const idx = getIJKComp(ijk, d);
  cell_t const dc = getDeltaCoord(d);
  geom_t clow = idx * dc;
  geom_t chigh = (idx + 1) * dc;
  return extents_t(clow, chigh);
  // return extents_t(m_bdys[cell-1],m_bdys[cell]);
}

template <typename cell_t, typename geometry_t, typename bdy_descriptor_t>
template <typename RNG_T>
typename Cartesian_3D<cell_t, geometry_t, bdy_descriptor_t>::Vector
Cartesian_3D<cell_t, geometry_t, bdy_descriptor_t>::sample_direction_isotropic(
    RNG_T & rng)
{
  geom_t const ctheta = geom_t(2) * rng.random() - geom_t(1);
  geom_t const phi = geom_t(2) * pi * rng.random();
  geom_t const stheta = std::sqrt(1 - ctheta * ctheta);
  geom_t const cphi = std::cos(phi);
  geom_t const sphi = std::sin(phi);
  Vector v;
  v.v = {{stheta * cphi, stheta * sphi, ctheta}};
  return v;
}

template <typename cell_t, typename geometry_t, typename bdy_descriptor_t>
template <typename RNG_T>
typename Cartesian_3D<cell_t, geometry_t, bdy_descriptor_t>::Vector
Cartesian_3D<cell_t, geometry_t, bdy_descriptor_t>::sample_position(
    RNG_T & rng,
    cell_t const cell) const
{
  Vector coord;
  extents_t xs = cell_extents(cell, X);
  extents_t ys = cell_extents(cell, Y);
  extents_t zs = cell_extents(cell, Z);
  coord.v[0] = xs.first + rng.random() * (xs.second - xs.first);

  coord.v[1] = ys.first + rng.random() * (ys.second - ys.first);

  coord.v[2] = zs.first + rng.random() * (zs.second - zs.first);
  return coord;
}  // sample_position

template <typename cell_t, typename geometry_t, typename bdy_descriptor_t>
cell_t
Cartesian_3D<cell_t, geometry_t, bdy_descriptor_t>::cell_across(
    cell_t const cell,
    face_t face,
    Vector const & /* p*/) const
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
    result = null_cell_;
  }
  return result;
}  // cell_across_face

template <typename cell_t, typename geometry_t, typename bdy_descriptor_t>
bool
Cartesian_3D<cell_t, geometry_t, bdy_descriptor_t>::isBoundary(
    cell_t const cidx,
    face_t const face) const
{
  bool isBound(false);
  ijk_T const ijk(mkIJK(cidx));
  cell_t const idx = getIJKComp(ijk, face);
  bool const ishigh(isHigh(face));
  if(0 == idx && !ishigh) { isBound = true; }
  else {
    cell_t ndir = getDimSize(face);
    if(ndir - 1 == idx && ishigh) isBound = true;
  }
  return isBound;
}

template <typename cell_t, typename geometry_t, typename bdy_descriptor_t>
Cartesian_3D<cell_t, geometry_t, bdy_descriptor_t>::Cartesian_3D(
    cell_t const nx,
    cell_t const ny,
    cell_t const nz,
    geom_t const dx,
    geom_t const dy,
    geom_t const dz,
    vec_bdy_descs const bds)
    : m_nx(nx),
      m_ny(ny),
      m_nz(nz),
      m_dx(dx),
      m_dy(dy),
      m_dz(dz),
      m_ncells(nx * ny * nz),
      m_bdy_descs(bds)
{
  nut::Require(m_bdy_descs.size() == 2 * nx * ny + 2 * nx * nz + 2 * ny * nz,
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
}  // ctor

template <typename cell_t, typename geometry_t, typename bdy_descriptor_t>
typename Cartesian_3D<cell_t, geometry_t, bdy_descriptor_t>::bdy_desc_t
Cartesian_3D<cell_t, geometry_t, bdy_descriptor_t>::getBC(
    face_t const face,
    cell_t const cidx) const
{
  cell_t const offset(getBCOffset(face, m_nx, m_ny, m_nz));
  ijk_T const ijk(mkIJK(cidx));

  dir_t dir = static_cast<dir_t>(face % 3);

  cell_t i1;
  cell_t i2;
  cell_t n1;
  switch(dir) {
    case X:
      i1 = ijk.j;
      i2 = ijk.k;
      n1 = m_ny;
      break;
    case Y:
      i1 = ijk.i;
      i2 = ijk.k;
      n1 = m_nx;
      break;
    case Z:
      i1 = ijk.i;
      i2 = ijk.j;
      n1 = m_nx;
      break;
  }
  cell_t const bcidx = getBCIndex(i1, i2, n1);
  cell_t const index = offset + bcidx;
  return m_bdy_descs[index];
}  // getBC

}  // namespace nut

// End of file
