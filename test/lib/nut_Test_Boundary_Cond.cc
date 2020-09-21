// nut_Test_Boundary_Cond.cc
// T. M. Kelley
// Nov 26, 2019
// (c) Copyright 2019 Triad National Security, all rights reserved

#include "Boundary_Cond.hh"
#include "Mesh_3D_Cartesian.hh"
#include "expect.hh"
#include "gtest/gtest.h"
#include "types.hh"  // boundary types

TEST(nut_Boundary_Condition, init)
{
  using Mesh_Iface_T = nut::Cartesian_Mesh_Interface;
  using Face = Mesh_Iface_T::face_handle_t;
  using BC_T = nut::Boundary_Cond<Face>;

  BC_T bcs;

  return;
}  //

TEST(nut_Boundary_Condition, 3D_examples)
{
  using Mesh_Iface_T = nut::Cartesian_Mesh_Interface;
  using Face = Mesh_Iface_T::face_handle_t;
  using BC_T = nut::Boundary_Cond<Face>;
  using desc_t = BC_T::desc_t;

  BC_T bcs;

  Face f0{0};
  EXPECT_FALSE(bcs.is_known_boundary(f0));
  bcs.set_boundary_type(f0, nut::bdy_types::VACUUM);
  EXPECT_TRUE(bcs.is_known_boundary(f0));
  desc_t d0 = bcs.get_boundary_type(f0);
  EXPECT_EQ(d0, nut::bdy_types::VACUUM);

  // Many faces should not be known to the boundary conditions
  for(uint32_t i = 1; i < 1000; ++i) {
    Face f{i};
    EXPECT_FALSE(bcs.is_known_boundary(f));
    bcs.set_boundary_type(f, nut::bdy_types::VACUUM);
    EXPECT_TRUE(bcs.is_known_boundary(f));
    desc_t d = bcs.get_boundary_type(f);
    EXPECT_EQ(d, nut::bdy_types::VACUUM);
  }
  return;
}  //

// End of file
