// Cartesian_Mesh.cc
// May 23, 2019
// (c) Copyright 2019 Triad National Security, all rights reserved

#include "Cartesian_Mesh.h"

namespace nut_mesh {

// initialize static data member
const Cartesian_Mesh::cell_handle_t Cartesian_Mesh::null_cell_{void_cell_idx};

const Cartesian_Mesh::face_handle_t Cartesian_Mesh::null_face_{void_cell_idx};
} // namespace nut_mesh

// End of file
