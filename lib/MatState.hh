// MatState.hh
// T. M. Kelley
// Jul 13, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

/* Intent is to aggregate object creation between fileio and app */

#ifndef MATSTATE_HH
#define MATSTATE_HH

#include "Assert.hh"
#include "Cell_Data.hh"
#include "Mesh.hh"
#include "constants.hh"
#include "fileio.hh"
#include <algorithm>

#ifdef HAVE_MURMELN
#include "murmeln/mesh_adaptors/Spherical_Mesh_Interface.h"
#endif

namespace nut {

// convert a row from the input file to a Cell_Data instance
template <typename fp_t, typename vector_t>
Cell_Data<fp_t, vector_t>
row_to_cell_data(MatStateRowP<fp_t, vector_t> const & row)
{
  fp_t const rho_p{(row.density / pmg)};
  fp_t const rho_e_minus{row.ye * rho_p};
  fp_t const rho_e_plus{0};
  fp_t const rho_A{0};
  fp_t const y_e{row.ye};
  fp_t const abar{0};

  fp_t const T_p{k_B * fp_t(1e9) * row.temperature};
  fp_t const T_e_minus{T_p};
  fp_t const T_e_plus{T_p};
  fp_t const l_nue{1e51 * (row.lnue_capture + row.lnue_pair)};
  fp_t const l_nueb{1e51 * (row.lnueb_capture + row.lnueb_pair)};
  fp_t const l_nux{1e51 * row.lnux_pair};
  vector_t v{row.velocity};

  return Cell_Data<fp_t, vector_t>{
      rho_p,     rho_e_minus, rho_e_plus, rho_A,  y_e,   abar, T_p,
      T_e_minus, T_e_plus,    l_nue,      l_nueb, l_nux, v};
};

// Convert all the rows to vector<Cell_Data>
template <typename fp_t, typename vector_t>
std::vector<Cell_Data<fp_t, vector_t>>
make_cell_data(std::vector<MatStateRowP<fp_t, vector_t>> const & rows)
{
  std::vector<Cell_Data<fp_t, vector_t>> cells{rows.size()};
  std::transform(rows.begin(), rows.end(), cells.begin(),
                 row_to_cell_data<fp_t, vector_t>);
  return cells;
};

/*!\brief: create a mesh within the given limits; identify the limiting
 * indices in the rows vector. The limiting indices use STL begin,end
 * convention. */
template <typename fp_t>
Spherical_1D_Mesh
rows_to_mesh(
    std::vector<MatStateRowP<fp_t, Spherical_1D_Mesh::Vector>> const & rows,
    fp_t const llimit,
    fp_t const ulimit,
    size_t & llimitIdx,
    size_t & ulimitIdx)
{
  using mesh_t = Spherical_1D_Mesh;
  using vector_t = mesh_t::Vector;
  using row_t = MatStateRowP<fp_t, vector_t>;
  // bool const  lims_ok = (ulimit > llimit);
  // Require(lims_ok , "rows_to_mesh: lower limit >= upper limit");
  size_t nrows = rows.size();
  std::vector<typename mesh_t::geom_t> bndsTmp(nrows + 1);
  size_t lIdx = 0;
  size_t uIdx = nrows;
  // get radii from rows
  std::vector<fp_t> rads(nrows);
  std::transform(rows.begin(), rows.end(), rads.begin(),
                 [](row_t const & row) { return row.radius; });
  bndsTmp[0] = rads[0] - (rads[1] - rads[0]) / 2;
  for(size_t i = 0; i < nrows - 1; ++i) {
    bndsTmp[i + 1] = rads[i] + (rads[i + 1] - rads[i]) / 2;
  }
  bndsTmp[nrows] = rads[nrows - 1] + (rads[nrows - 1] - rads[nrows - 2]) / 2;

  // would use STL algs here, but need a actual indices
  // find limiting indices: want the greatest index such that
  // bndsTmp[i] < llimit
  while(bndsTmp[lIdx + 1] < llimit && lIdx < nrows) { lIdx++; }
  if(lIdx == nrows) {
    std::stringstream errstr;
    errstr << "rows_to_mesh: lower bound (" << llimit
           << ") not within range of input (max input radius: "
           << rads[nrows - 1] << ", which comes out to upper bdy of "
           << bndsTmp[nrows] << std::endl;
    throw(std::runtime_error(errstr.str()));
  }
  // now look for the least upper index such that
  // bndsTmp[uIdx] > ulimit
  uIdx = lIdx;
  while(bndsTmp[uIdx] < ulimit && uIdx < nrows) { uIdx++; }

  // outputs
  ulimitIdx = uIdx;
  llimitIdx = lIdx;
  size_t const ncells = uIdx - lIdx;
  std::vector<typename mesh_t::geom_t> bounds(ncells + 1);
  std::vector<typename mesh_t::bdy_desc_t> descs(ncells + 1);

  std::copy(&bndsTmp[lIdx], &bndsTmp[uIdx + 1], bounds.begin());
  descs[0] = bdy_types::REFLECTIVE;
  descs[ncells] = bdy_types::VACUUM;
  for(size_t i = 1; i < ncells; ++i) { descs[i] = bdy_types::NONE; }
  return mesh_t(bounds, descs);
}  // rows_to_mesh

#ifdef HAVE_MURMELN
/*!\brief: create a mesh within the given limits; identify the limiting
 * indices in the rows vector. The limiting indices use STL begin,end
 * convention. */
template <typename fp_t>
murmeln_mesh::Spherical_1D_Mesh
rows_to_murmeln_mesh(
    std::vector<
        MatStateRowP<fp_t, murmeln::Spherical_Mesh_Interface::Vector>> const &
        rows,
    fp_t const llimit,
    fp_t const ulimit,
    size_t & llimitIdx,
    size_t & ulimitIdx)
{
  using murmeln_mesh::Spherical_1D_Mesh;
  using mesh_t = Spherical_1D_Mesh;
  using vector_t = mesh_t::Vector;
  using row_t = MatStateRowP<fp_t, vector_t>;
  // bool const  lims_ok = (ulimit > llimit);
  // Require(lims_ok , "rows_to_mesh: lower limit >= upper limit");
  size_t nrows = rows.size();
  std::vector<fp_t> bndsTmp(nrows + 1);
  size_t lIdx = 0;
  size_t uIdx = nrows;
  // get radii from rows
  std::vector<Spherical_1D_Mesh::Geom_T> rads(nrows);
  std::transform(rows.begin(), rows.end(), rads.begin(),
                 [](row_t const & row) { return row.radius; });
  bndsTmp[0] = rads[0] - (rads[1] - rads[0]) / 2;
  for(size_t i = 0; i < nrows - 1; ++i) {
    bndsTmp[i + 1] = rads[i] + (rads[i + 1] - rads[i]) / 2;
  }
  bndsTmp[nrows] = rads[nrows - 1] + (rads[nrows - 1] - rads[nrows - 2]) / 2;

  // would use STL algs here, but need a actual indices
  // find limiting indices: want the greatest index such that
  // bndsTmp[i] < llimit
  while(bndsTmp[lIdx + 1] < llimit && lIdx < nrows) { lIdx++; }
  if(lIdx == nrows) {
    std::stringstream errstr;
    errstr << "rows_to_mesh: lower bound (" << llimit
           << ") not within range of input (max input radius: "
           << rads[nrows - 1] << ", which comes out to upper bdy of "
           << bndsTmp[nrows] << std::endl;
    throw(std::runtime_error(errstr.str()));
  }
  // now look for the least upper index such that
  // bndsTmp[uIdx] > ulimit
  uIdx = lIdx;
  while(bndsTmp[uIdx] < ulimit && uIdx < nrows) { uIdx++; }

  // outputs
  ulimitIdx = uIdx;
  llimitIdx = lIdx;
  size_t const ncells = uIdx - lIdx;
  std::vector<Spherical_1D_Mesh::Geom_T> bounds(ncells + 1);

  std::copy(&bndsTmp[lIdx], &bndsTmp[uIdx + 1], bounds.begin());
  return mesh_t(std::move(bounds));
}  // rows_to_murmeln_mesh
#endif
// HAVE_MURMELN

}  // namespace nut

#endif  // include guard

// End of file
