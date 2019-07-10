// Mesh_3D_Cartesian.hh
// T. M. Kelley
// May 07, 2019
// (c) Copyright 2019 Triad National Security, all rights reserved

#pragma once

#ifdef HAVE_MURMELN
#include "murmeln/mesh_adaptors/Cartesian_Mesh_Interface.h"

namespace nut {

inline murmeln::Cartesian_Mesh_Interface
make_cartesian_mesh_interface()
{
  using Index_T = murmeln::Cartesian_Mesh_Interface::Index_T;
  using Geom_T = murmeln::Cartesian_Mesh_Interface::Geom_T;
  Index_T const nx = 3, ny = 5, nz = 7;
  Geom_T const dx = 1.0, dy = 2.0, dz = 3.0;
  // vacuum bdy conds are default
  murmeln::Cartesian_Mesh_Interface mesh(nx, ny, nz, dx, dy, dz);
  return mesh;
}

}  // namespace nut

#endif
// #ifdef HAVE_MURMELN

// End of file
