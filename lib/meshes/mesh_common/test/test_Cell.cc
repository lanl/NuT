// test_cell.cc
// Jan 10, 2019
// (c) Copyright 2019 LANSLLC, all rights reserved

#include "mesh_common/Cell.h"
#include "gtest/gtest.h"

TEST(default_mesh_cell, instantiate) {
  using namespace nut_mesh;
  Cell c(420u);

  EXPECT_EQ(c.as_id(), 420u);

} // TEST(default_mesh,cell){

// End of file
