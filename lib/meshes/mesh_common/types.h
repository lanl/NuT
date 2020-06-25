// default_mesh_types.h
// Dec 21, 2018
// (c) Copyright 2018 LANSLLC, all rights reserved

#pragma once

#include <cstdint>
#include <limits>

namespace murmeln_mesh {

using index_t = ::uint64_t;

index_t constexpr max_index_t = std::numeric_limits<index_t>::max();

using cell_index_t = index_t;

using geom_t = double;

} // namespace murmeln_mesh

// End of file
