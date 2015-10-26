// mpi.hh
// T. M. Kelley
// Jun 20, 2012
// (c) Copyright 2012 LANSLLC, all rights reserved


#ifndef MPI_HELPERS_HH
#define MPI_HELPERS_HH

#include "Assert.hh"
#include <stdint.h>

/** types and helpers for MPI */
namespace nut
{

    struct Rank
    {
        Rank(uint32_t rank) : r(rank) {}
        uint32_t fromRank() const {return r;}
        uint32_t r;
    };

    struct CommSz
    {
        CommSz(uint32_t commSz) : s(commSz) {
            dbc::GreaterThan(s,0u,"communicator size");
        }
        uint32_t fromCommSz() const {return s;}
        uint32_t s;
    };

} // nut::


#endif // include guard


// End of file
