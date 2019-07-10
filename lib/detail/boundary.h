// boundary.h
// Jun 25, 2019
// (c) Copyright 2019 Triad National Security, all rights reserved

#pragma once

#include "types.hh"

namespace nut::boundary_types {

// Cell boundary maps
enum BDY_TYPE : index_t {
  NONE = 0,
  CELL = NONE,
  VACUUM = 1,
  REFLECTIVE = 2,
  PERIODIC = 3,
  PROCESSOR = 4
};

inline std::string
to_string(BDY_TYPE b) {
  std::string s;
  switch (b) {
    case NONE:
      s = "NONE";
      break;
    case VACUUM:
      s = "VACUUM";
      break;
    case REFLECTIVE:
      s = "REFLECTIVE";
      break;
    case PERIODIC:
      s = "PERIODIC";
      break;
    case PROCESSOR:
      s = "PROCESSOR";
      break;
  }
  return s;
}

}  // namespace nut::boundary

// End of file
