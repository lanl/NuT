// Mesh_3D_Cartesian.hh
// T. M. Kelley
// May 07, 2019
// (c) Copyright 2019 Triad National Security, all rights reserved

#pragma once

#ifdef HAVE_MURMELN

#include "murmeln/mesh_adaptors/Cartesian_Mesh_Interface.h"

namespace nut {

inline murmeln_mesh::Cartesian_Mesh
make_cartesian_mesh()
{
  using Index_T = murmeln::Cartesian_Mesh_Interface::Index_T;
  using Geom_T = murmeln::Cartesian_Mesh_Interface::Geom_T;
  Index_T const nx = 3;
  Index_T const ny = 5;
  Index_T const nz = 7;
  Geom_T const dx = 1.0;
  Geom_T const dy = 2.0;
  Geom_T const dz = 3.0;
  // vacuum bdy conds are default
  murmeln_mesh::Cartesian_Mesh mesh(nx, ny, nz, dx, dy, dz);
  return mesh;
}

}  // namespace nut

#else

#include "Cartesian_Mesh.h"

namespace nut {
inline nut::Cartesian_Mesh
make_cartesian_mesh()
{
  using Index_T = Cartesian_Mesh::Index_T;
  using Geom_T = Cartesian_Mesh::Geom_T;
  Index_T const nx = 3;
  Index_T const ny = 5;
  Index_T const nz = 7;
  Geom_T const dx = 1.0;
  Geom_T const dy = 2.0;
  Geom_T const dz = 3.0;
  // vacuum bdy conds are default
  Cartesian_Mesh mesh(nx, ny, nz, dx, dy, dz);
  return mesh;
}  // make_cartesian_mesh

}  // namespace nut

#endif
// #ifdef HAVE_MURMELN

// End of file
