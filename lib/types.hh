// types.hh
// T. M. Kelley
// Jan 03, 2011
// Header for types
// (c) Copyright 2011 LANSLLC all rights reserved.

#pragma once

#include "Mesh_3D_Cartesian.hh"  // TO DO this is very likely circular.
                                 // likely solution: include murmeln header directly
#include "RNG.hh"
#include "Vec3D.hh"
#include <iterator>
#include <stdexcept>
#include <string>
#include <type_traits>  // std::is_integral
#include <utility>      // std::pair
#include <vector>
#include <stdint.h>

namespace nut {

namespace bdy_types {

#ifdef HAVE_MURMELN
using namespace murmeln_mesh::boundary;
using descriptor = BDY_TYPE;
// using descriptor = murmeln::Cartesian_Mesh_Interface::Boundary;
// enum descriptor {
//   NONE = murmeln::Cartesian_Mesh_Interface::Boundary::NONE,
//   CELL = murmeln::Cartesian_Mesh_Interface::Boundary::CELL,
//   VACUUM = murmeln::Cartesian_Mesh_Interface::Boundary::VACUUM,
//   REFLECTIVE = murmeln::Cartesian_Mesh_Interface::Boundary::REFLECTIVE,
//   PERIODIC = murmeln::Cartesian_Mesh_Interface::Boundary::PERIODIC,
//   PROCESSOR = murmeln::Cartesian_Mesh_Interface::Boundary::PROCESSOR,
// };
#else
enum descriptor {
  NONE = 0,
  CELL = NONE,
  VACUUM = 1,
  REFLECTIVE = 2,
  PERIODIC = 3,
  PROCESSOR = 4,
};
#endif  // HAVE_MURMELN

}  // namespace bdy_types

using geom_t = double; /*^ type for geometry calculations */

/** \brief Our type for geometric vectors */
template <size_t dim>
using vec_t = Vec_T<geom_t, dim>;

using vg = std::vector<geom_t>;

using cell_t = uint32_t; /*^ cell index. Use cell_t(-1) for null cell. */

using id_t = uint32_t; /*^ particle index */

using vid = std::vector<id_t>;

using group_t = uint32_t;

using cntr_t = uint32_t;

using arg_error = std::invalid_argument;

using sz_o_it = std::ostream_iterator<size_t>;

using ge_o_it = std::ostream_iterator<geom_t>;

using rng_t = nut::Philox4x32_RNG;

using ctr_t = rng_t::ctr_t;

using key_t = rng_t::key_t;

using seed_t = rng_t::seed_t;

/* Neutrino species */
enum Species {
  nu_e,
  nu_e_bar,
  nu_x,
  nu_mu,
  nu_mu_bar,
  nu_tau,
  nu_tau_bar
};  // Species;

std::string
species_name(Species const s);

Species
species_type(std::string const & s);

seed_t
species_seed(Species const s);

}  // namespace nut

// End of file
