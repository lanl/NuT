// Particle.hh
// T. M. Kelley
// Dec 21, 2010
// Header for Particle
// (c) Copyright 2010 LANSLLC all rights reserved.

#ifndef PARTICLE_H
#define PARTICLE_H

#include "Fates.hh"
#include "types.hh"
#include "utilities.hh"
#include <algorithm>
#include <ostream>

namespace nut {

/**\class Particle
 * \tparam fpt: floating point type
 * \tparam rngt: random number generator type
 * \tparam vectort: space vector type
 * \tparam cellt: cell index type
 */
template <typename fpt, typename rngt, typename vectort>
struct Particle {
  // types and class constants
  // static const uint32_t dim = Dim;
  using fp_t = fpt;
  using rng_t = rngt;
  using vec_t = vectort;
  using cell_t = uint32_t;
  static constexpr size_t dim = vectort::dim;

  vec_t x;
  vec_t omega;
  fp_t e;
  fp_t t;
  fp_t weight;
  cell_t cell;
  rng_t rng;
  Species species;
  bool alive;
  Fates fate;

  // ctor
  Particle(vec_t x_,
           vec_t omega_,
           fp_t e_,
           fp_t t_,
           fp_t weight_,
           cell_t cell_,
           rng_t rng_,
           Species s_)
      : x(x_),
        omega(omega_),
        e(e_),
        t(t_),
        weight(weight_),
        cell(cell_),
        rng(rng_),
        species(s_),
        alive(true),
        fate(NOT_DEAD_YET)
  {
  }

  Particle() {}

  void kill() { alive = false; }

  bool operator==(Particle const & p) const
  {
    bool passed = arrays_eq(x, p.x) && omega == p.omega && e == p.e &&
                  t == p.t && weight == p.weight && cell == p.cell &&
                  rng == p.rng && species == p.species && alive == p.alive;
    return passed;
  }

};  // Particle

}  // namespace nut

#endif

// version
// $Id$

// End of file
