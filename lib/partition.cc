// partition.cc
// T. M. Kelley
// Jun 19, 2012
// (c) Copyright 2012 LANSLLC, all rights reserved

#include "partition.hh"
#include <iostream>

namespace nut {

void
getChunks(uint32_t const rank,
          uint32_t const commSz,
          uint32_t const nc,
          ChunkIdVec & chunks)
{
  uint32_t nr = (nc / commSz) + ((rank < nc % commSz) ? 1 : 0);
  std::cout << "getChunks: n_chks = " << nc << ", commSz = " << commSz
            << ", rank = " << rank << ", nr = " << nr
            << ", nc / commSz = " << (nc / commSz) << std::endl;
  // uint32_t n = nc / commSz;
  chunks.resize(nr);
  if(nr > 0) {
    for(uint32_t i = 0; i < nr; ++i) { chunks[i] = rank + i * commSz; }
  }
  return;
}  // getChunks

}  // namespace nut

// End of file
