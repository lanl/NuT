// T. M. Kelley (c) 2011 LANS LLC

#include "Particle.hh"
#include "RNG.hh"
#include "gtest/gtest.h"
#include "meshes/mesh_common/Vector.h"

TEST(nut_Particle, instantiation)
{
  // for a very simple test, we can use a totally bogus RNG
  typedef uint32_t rng_t;
  typedef double fp_t;

  using Vector = murmeln_mesh::Vector;

  typedef nut::Particle<fp_t, rng_t, Vector> part_t;

  fp_t x(1.0);
  fp_t omega(1.0);
  fp_t e(1.0);
  fp_t t(1.0);
  fp_t wt(1.0);
  nut::cell_t cell(1);
  nut::Species s(nut::nu_e);

  rng_t rng(42);

  part_t particle({x}, {omega}, e, t, wt, cell, rng, s);
  EXPECT_TRUE(true);
  return;
}

// End of file
