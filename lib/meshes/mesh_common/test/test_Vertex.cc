// test_Vertex.cc
// Jan 10, 2019
// (c) Copyright 2019 LANSLLC, all rights reserved

#include "mesh_common/Vertex.h"
#include "gtest/gtest.h"

TEST(default_mesh_Vertex, instantiate) {
  using namespace murmeln_mesh;
  Vertex v(67);
  EXPECT_EQ(v.as_id(), 67);
}

// End of file
