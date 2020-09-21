// test_Face.cc
// Jan 10, 2019
// (c) Copyright 2019 LANSLLC, all rights reserved

#include "mesh_common/Face.h"
#include "gtest/gtest.h"

TEST(default_mesh_Face, instantiate) {
  using namespace nut_mesh;
  Face f(56);
  EXPECT_EQ(f.as_id(), 56);
}

TEST(default_mesh_Face, hashable) {
  using namespace nut_mesh;
  {
    Spherical_1D_Face f(56);
    std::hash<Spherical_1D_Face> fhash{};
    size_t hashf{fhash(f)};
    std::hash<index_t> ihash{};
    size_t hashi{ihash(56)};
    EXPECT_EQ(hashi, hashf);
  }
  {
    Cartesian_Face f(56);
    std::hash<Cartesian_Face> fhash{};
    size_t hashf{fhash(f)};
    std::hash<index_t> ihash{};
    size_t hashi{ihash(56)};
    EXPECT_EQ(hashi, hashf);
  }
  return;
}

// End of file
