// test_Edge.cc
// Jan 10, 2019
// (c) Copyright 2019 LANSLLC, all rights reserved

#include "mesh_common/Edge.h"
#include "gtest/gtest.h"

TEST(default_mesh_Edge, instantiate) {
  using namespace nut_mesh;
  Edge e(45);
  EXPECT_EQ(e.as_id(), 45);
}

// End of file
