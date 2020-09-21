/**
 * \file Header file of vertex class
 */

#pragma once

#include "mesh_common/Mesh_Element.h"

namespace nut_mesh {

/**
 * @brief Geometry vertex class
 *
 * Vertex class, containing all information of geometry vertices.
 */
struct Vertex : public Mesh_Element {
public:
  explicit Vertex(index_t id) : Mesh_Element(id) {}

}; // Vertex

} // namespace nut_mesh
