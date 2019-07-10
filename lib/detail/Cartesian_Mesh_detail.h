// Cartesian_Mesh_detail.h
// Apr 30, 2019
// (c) Copyright 2019 Triad National Security, all rights reserved
/* Free functions that help implement the Cartesian mesh */

#pragma once

#include "types.hh"
#include <tuple>

namespace nut {

using CIndices = std::tuple<index_t, index_t, index_t>;

// Free functions that help define Cartesian relationships

/**\brief Convert a Cartesian triple of indices to a unique linear index.
 * \param ix, iy, iz: Cartesian indices of cell
 * \param nx, ny: Number of cells in X- and Y-dimension
 */
inline index_t
cartesian_to_linear(index_t ix, index_t iy, index_t iz, index_t nx, index_t ny)
{
  return ix + (iy + iz * ny) * nx;
}

/**\brief Convert a linear 3D index into a Cartesian triple of indices. */
inline CIndices
linear_to_cartesian(index_t const idx, index_t nx, index_t ny)
{
  index_t const ix = idx % nx;
  index_t const rem_xy = idx / nx;
  index_t const iy = rem_xy % ny;
  index_t const iz = rem_xy / ny;
  return {ix, iy, iz};
}

/**\brief The number of faces perpendicular to the X direction/
 *\param nx: number of x cells in mesh
 *\param ny: number of y cells in mesh
 *\param nz: number of z cells in mesh
 */
inline index_t
num_yz_faces(index_t nx, index_t ny, index_t nz)
{
  return ny * nz * (nx + 1);
}

/**\brief The number of faces perpendicular to the Y direction
 *\param nx: number of x cells in mesh
 *\param ny: number of y cells in mesh
 *\param nz: number of z cells in mesh
 */
inline index_t
num_xz_faces(index_t nx, index_t ny, index_t nz)
{
  return nx * nz * (ny + 1);
}

/**\brief The number of faces perpendicular to the Z direction
 *\param nx: number of x cells in mesh
 *\param ny: number of y cells in mesh
 *\param nz: number of z cells in mesh
 */
inline index_t
num_xy_faces(index_t nx, index_t ny, index_t nz)
{
  return nx * ny * (nz + 1);
}

}  // namespace nut

// End of file
