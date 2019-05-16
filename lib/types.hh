// types.hh
// T. M. Kelley
// Jan 03, 2011
// Header for types
// (c) Copyright 2011 LANSLLC all rights reserved.

#ifndef TYPES_H
#define TYPES_H

#include "RNG.hh"
#include "Vec3D.hh"
#include <iterator>
#include <stdexcept>
#include <string>
#include <utility>  // std::pair
#include <vector>
#include <stdint.h>

namespace nut {
namespace bdy_types {
enum descriptor {
  V, /** Vacuum */
  R, /** Reflective */
  T  /** Transmissive */
};
}

namespace events {
enum Event {
  // physics events
  collision = 0,
  nucleon_abs = 1,
  nucleon_elastic_scatter = 2,
  electron_scatter = 3,
  positron_scatter = 4,
  nu_e_annhilation = 5,
  nu_x_annhilation = 6,
  // # compute events
  boundary = 100,
  cell_low_x_boundary = 101,
  cell_high_x_boundary = 102,
  escape = 110,
  reflect = 111,
  step_end = 115,
  weight_cutoff = 120,
  // testing only
  null = 9999
};
}  // namespace events

using geom_t = double; /*^ type for geometry calculations */

/** \brief Our type for geometric vectors */
template <size_t dim>
using vec_t = Vec_T<geom_t, dim>;

using vg = std::vector<geom_t>;

using event_n_dist = std::pair<events::Event, geom_t>;

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

namespace events {
std::string
event_name(Event const & e);
}  // namespace events

}  // namespace nut

#endif

// version
// $Id$

// End of file
