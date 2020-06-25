// test_Mesh_Element.cc
// Jan 10, 2019
// (c) Copyright 2019 LANSLLC, all rights reserved

#include "mesh_common/Mesh_Element.h"
#include "gtest/gtest.h"

TEST(default_mesh, Mesh_Element) {
  using murmeln_mesh::Mesh_Element;
  Mesh_Element e(42u);
  EXPECT_EQ(e.as_id(), 42u);

  Mesh_Element e1(43u);
  EXPECT_FALSE(e1 == e);
  EXPECT_TRUE(e == e);
  EXPECT_TRUE(e1 == e1);
  EXPECT_TRUE(e1 != e);
}

// End of file
