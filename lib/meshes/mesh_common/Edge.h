/**
 * \file for the Edge class
 */

#pragma once

#include <vector>

#include "mesh_common/Mesh_Element.h"

namespace murmeln_mesh {

/**
 * @class An edge
 */
struct Edge : public Mesh_Element {

  explicit Edge(index_t id) : Mesh_Element(id) {}
};

} // namespace murmeln_mesh
