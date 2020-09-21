// Spherical_Mesh.cc
// Jun 27, 2019
// (c) Copyright 2019 Triad National Security, all rights reserved

#include "spherical/Spherical_Mesh.h"
namespace nut_mesh {

/* Note, using this value to distinguish from accidental roll-over cases. */
const Spherical_1D_Mesh::Face Spherical_1D_Mesh::null_face_{0xFFFF'FFFC};

const Spherical_1D_Mesh::Cell Spherical_1D_Mesh::null_cell_{0xFFFF'FFFC};
} // namespace nut_mesh

// End of file
